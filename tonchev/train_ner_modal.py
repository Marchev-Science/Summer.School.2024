from pathlib import Path
from modal import Image, App, Volume
import modal
from datasets import load_dataset
from transformers import (AutoTokenizer,
                          DataCollatorForTokenClassification,
                          AutoModelForTokenClassification,
                          TrainingArguments,
                          Trainer,
                          )
import numpy as np
import evaluate
import warnings
warnings.filterwarnings("ignore")

NER_TAGS_COL: str = 'tags'
OUTPUT_DIR: str = '/roberta-ner'

# modal stuff
MODEL_DIR = "/model"


debian_cuda_image = (
    Image.debian_slim(python_version='3.11')
    .pip_install("torch", gpu='any')
    .pip_install('transformers',
                 'datasets',
                 'evaluate',
                 'accelerate',
                 'seqeval',
                 'nvidia-ml-py3',
                 )
)

# cuda_torch_image = (
#     Image.micromamba(python_version='3.11')
#     .micromamba_install('pytorch', 'pytorch-cuda=12.1', gpu='any', channels=['pytorch', 'nvidia', 'conda-forge'])
#     .pip_install('transformers',
#                  'accelerate',
#                  'datasets',
#                  'evaluate',
#                  'seqeval',
#                  'nvidia-ml-py3',
#                  )
# )

volume = Volume.from_name("model-volume", create_if_missing=True)

app = App("ner-train")


def print_gpu_utilization():
    from pynvml import nvmlInit, nvmlDeviceGetHandleByIndex, nvmlDeviceGetMemoryInfo
    nvmlInit()
    handle = nvmlDeviceGetHandleByIndex(0)
    info = nvmlDeviceGetMemoryInfo(handle)
    print(f"GPU memory occupied: {info.used//1024**2} MB.")


def print_summary(result):
    print(f"Time: {result.metrics['train_runtime']:.2f}")
    print(f"Samples/second: {result.metrics['train_samples_per_second']:.2f}")
    print_gpu_utilization()

def info():
    if not modal.is_local():
        import torch
        print("I am on modal")
        print(list(Path('.').iterdir()))
        print(Path('.').absolute())
        print(f'torch device: {torch.cuda.get_device_name()}')
        print(f'torch version: {torch.__version__}')
        print(f'cuda version: {torch.version.cuda}')
        print_gpu_utilization()


# ner utils
def relable(examples, original_to_ours):
    """Relables the ner tags to more convenient order"""
    new_tags = []
    for tags in examples[NER_TAGS_COL]:
        new_tags.append([original_to_ours[item] for item in tags])
    examples[NER_TAGS_COL] = new_tags
    return examples


def align_labels_with_tokens(labels, word_ids):
    new_labels = []
    current_word = None
    for word_id in word_ids:
        if word_id != current_word:
            # Start of a new word!
            current_word = word_id
            label = -100 if word_id is None else labels[word_id]
            new_labels.append(label)
        elif word_id is None:
            # Special token
            new_labels.append(-100)
        else:
            # Same word as previous token
            label = labels[word_id]
            # If the label is B-XXX we change it to I-XXX
            if label % 2 == 1:
                label += 1
            new_labels.append(label)

    return new_labels


def tokenize_and_align_labels(examples, tokenizer, is_split_into_words=True):
    tokenized_inputs = tokenizer(
        examples["tokens"], truncation=True, is_split_into_words=is_split_into_words
    )
    all_labels = examples[NER_TAGS_COL]
    new_labels = []
    for i, labels in enumerate(all_labels):
        word_ids = tokenized_inputs.word_ids(i)
        new_labels.append(align_labels_with_tokens(labels, word_ids))

    tokenized_inputs["labels"] = new_labels
    return tokenized_inputs

def ner_metrics_factory(id2label: dict | list, module: str = "seqeval"):

    metric = evaluate.load(module)

    def compute_metrics(eval_preds):

        print_gpu_utilization()

        logits, labels = eval_preds
        predictions = np.argmax(logits, axis=-1)

        # Remove ignored index (special tokens) and convert to labels
        true_labels = [[id2label[l] for l in label if l != -100] for label in labels]
        true_predictions = [
            [id2label[p] for (p, l) in zip(prediction, label) if l != -100]
            for prediction, label in zip(predictions, labels)
        ]
        all_metrics = metric.compute(predictions=true_predictions, references=true_labels)
        return {
            "precision": all_metrics["overall_precision"],
            "recall": all_metrics["overall_recall"],
            "f1": all_metrics["overall_f1"],
            "accuracy": all_metrics["overall_accuracy"],
        }

    return compute_metrics

@app.function(image=debian_cuda_image,
              gpu='a10g',
              timeout=900,
              volumes={MODEL_DIR: volume}
              )
def train():

    info()

    raw_datasets = load_dataset("tner/mit_restaurant")
    model_checkpoint = "FacebookAI/roberta-base"

    original_label_to_id = {
        "O": 0,
        "B-Rating": 1,
        "I-Rating": 2,
        "B-Amenity": 3,
        "I-Amenity": 4,
        "B-Location": 5,
        "I-Location": 6,
        "B-Restaurant_Name": 7,
        "I-Restaurant_Name": 8,
        "B-Price": 9,
        "B-Hours": 10,
        "I-Hours": 11,
        "B-Dish": 12,
        "I-Dish": 13,
        "B-Cuisine": 14,
        "I-Price": 15,
        "I-Cuisine": 16
    }

    label2id = {
        "O": 0,
        "B-Rating": 1,
        "I-Rating": 2,
        "B-Amenity": 3,
        "I-Amenity": 4,
        "B-Location": 5,
        "I-Location": 6,
        "B-Restaurant_Name": 7,
        "I-Restaurant_Name": 8,
        "B-Price": 9,
        "I-Price": 10,
        "B-Hours": 11,
        "I-Hours": 12,
        "B-Dish": 13,
        "I-Dish": 14,
        "B-Cuisine": 15,
        "I-Cuisine": 16
    }
    id2label = {i: name for name, i in label2id.items()}

    original_to_ours = {original_label_to_id[k]: label2id[k] for k in label2id}

    tokenizer = AutoTokenizer.from_pretrained(model_checkpoint, add_prefix_space=True)
    model = AutoModelForTokenClassification.from_pretrained(
        model_checkpoint,
        id2label=id2label,
        label2id=label2id,
    )
    print('model is init')
    # model.to('cuda')
    print_gpu_utilization()

    # import sys
    # sys.exit()
    # prepare datasets
    raw_datasets = raw_datasets.map(relable, batched=True, fn_kwargs={'original_to_ours': original_to_ours})
    tokenized_datasets = raw_datasets.map(tokenize_and_align_labels,
                                          fn_kwargs={'tokenizer': tokenizer},
                                          batched=True,
                                          remove_columns=raw_datasets['train'].column_names)

    data_collator = DataCollatorForTokenClassification(tokenizer=tokenizer)

    # create the compute metric fn
    compute_metrics = ner_metrics_factory(id2label)

    # set training args
    args = TrainingArguments(
        MODEL_DIR + OUTPUT_DIR,
        eval_strategy='epoch',
        save_strategy='no',   # 'no', #epoch
        learning_rate=2e-5,
        num_train_epochs=5,
        weight_decay=0.01,
        push_to_hub=False,
        use_cpu=False,
        per_device_train_batch_size=300,
        logging_strategy='epoch',
    )

    trainer = Trainer(
        model=model,
        args=args,
        train_dataset=tokenized_datasets["train"],    # noqa .select(range(100)),
        eval_dataset=tokenized_datasets["validation"],   # noqa .select(range(100)),
        data_collator=data_collator,
        compute_metrics=compute_metrics,
        tokenizer=tokenizer,
    )

    trainer.train()
    result = trainer.train()
    print_summary(result)

    trainer.save_model(MODEL_DIR + OUTPUT_DIR)

@app.local_entrypoint()
def main():

    # run the function remotely on Modal
    train.remote()

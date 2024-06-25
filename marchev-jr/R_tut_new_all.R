# Учебна Сесия 1: Въведение и Основен Синтаксис в R

# Инсталиране
# https://cloud.r-project.org/
# https://posit.co/download/rstudio-desktop/
#
# https://cran.r-project.org/bin/windows/base/R-4.3.2-win.exe
# https://download1.rstudio.org/electron/windows/RStudio-2023.12.0-369.exe

# Интерфейс в R Studio
#
# 1. **Конзола (Console)**
# - Намира се в ляво долу (в основния изглед).
# - Тук се изпълнява R кодът в реално време.
# - Можете да въвеждате команди директно в конзолата и да виждате резултатите от тяхното изпълнение.
# 
# 2. **Скрипт Редактор (Script Editor)**
# - Разположен в горната лява част.
# - Тук можете да пишете и редактирате R скриптове (файлове с разширение `.R`).
# - Поддържа функции като синтактично оцветяване, автоматично подравняване и маркиране на скоби.
#
# 3. Работно пространство (Environment) на RStudio:
# - В дясната горна част на RStudio ще намерите раздела "Работно пространство" (Environment). 
# - Можете да видите текущите обекти и променливи във вашия проект.
# - В този раздел се извеждат променливите и техните стойности
# - Изчистване на всички обекти и променливи
# rm(list = ls())
# gc()
#
# 4. **Прозорец  и Файлове (History and Files)**
# - Разположен в долната дясна част.
# - Прозорецът за файлове (Files) позволява достъп до файловата система, показва файлове и папки.
# - Сетване на работната директория
# setwd("c:/my documents/work")
# setwd("c:\\my documents\\work")
# 
# 5. **Прозорец за История, Пакети, Помощ и Графики (Packages, Help, Plots)**
# - Намира се в долната дясна част.
# - Историята (History) показва списък на предишни команди.
# - Пакетите (Packages) показват инсталираните R пакети и дават възможност за инсталиране/деинсталиране.
# - Помощта (Help) предоставя документация за функциите на R.
# - Графиките (Plots) показват визуализации, създадени от R кода.
#
# 6. Ето някои от най-често използваните бързи клавиши в RStudio:
# - За да изпълните текущия ред или подсекция в скрипта, използвайте `Ctrl + Enter`.
# - За да коментирате или разкоментирате избрания ред, използвайте `Ctrl + Shift + C`.
# - За да запазите текущо редактирания файл, използвайте `Ctrl + S`.

# Операция 1: Декларация на променливи
# Пример 1: Деклариране на числова променлива
x <- 5
# Пример 2: Деклариране на символен низ
y <- "Здравей, R!"
# Пример 3: Деклариране на булева променлива
z <- TRUE

# Операция 2: Основни математически операции
# Пример 1: Събиране
a <- 3 + 7
# Пример 2: Умножение
b <- 8 * 2
# Пример 3: Деление
c <- 15 / 3

# Операция 3: Използване на логически оператори
# Пример 1: Логическо И (AND)
and_result <- (5 > 3) & (3 < 4)
# Пример 2: Логическо ИЛИ (OR)
or_result <- (5 < 3) | (3 == 3)
# Пример 3: Отрицание (NOT)
not_result <- !(5 == 5)

# Операция 4: Създаване и използване на вектори
# Пример 1: Създаване на числов вектор
vec_num <- c(1, 2, 3)
# Пример 2: Създаване на символен вектор
vec_char <- c("ябълка", "банан", "портокал")
# Пример 3: Създаване на булев вектор
vec_bool <- c(TRUE, FALSE, TRUE)

# Операция 5: Индексиране и извличане на елементи от вектор
# Пример 1: Извличане на първия елемент от вектор
first_elem <- vec_num[1]
# Пример 2: Извличане на последния елемент от вектор
last_elem <- vec_char[length(vec_char)]
# Пример 3: Извличане на всички елементи, освен първия
all_but_first <- vec_num[-1]

# Операция 6: Използване на матрици
# Пример 1: Създаване на матрица
matrix_example <- matrix(1:9, nrow = 3)
# Пример 2: Достъп до елементи на матрицата
element <- matrix_example[2, 3]
# Пример 3: Промяна на структурата на матрицата
reshaped_matrix <- matrix(matrix_example, nrow = 1)
# Пример 4: Транспониране на матрица
transposed_matrix <- t(matrix)
# Пример 5: Умножение на матрици
result_matrix <- matrix %*% t(matrix)
# Пример 6: Умножение на матрица с число
scalar_multiply <- matrix1 * 2

# Операция 7: Използване на условни оператори
# Пример 1: Проверка дали число е по-голямо от друго
is_greater <- 5 > 3
# Пример 2: Проверка за равенство между две числа
is_equal <- (2 + 2) == 4
# Пример 3: Проверка дали елемент съществува във вектор
is_in_vector <- "ябълка" %in% vec_char # грешен отговор
is_in_vector <- "ябълка" %in% unlist(vec_char)

# Операция 8: Използване на функции за агрегация
# Пример 1: Намиране на средно аритметично на вектор
mean_num <- mean(vec_num)
# Пример 2: Намиране на максимална стойност във вектор
max_num <- max(vec_num)
# Пример 3: Намиране на минимална стойност във вектор
min_num <- min(vec_num)

# Операция 9: Използване на string функции
# Пример 1: Обединение на стрингове
full_name <- paste("Иван", "Иванов")
# Пример 2: Разделяне на стрингове
name_parts <- strsplit(full_name, " ")[[1]]
# Пример 3: Замяна в стрингове
greet <- gsub("Иван", "Петър", "Здравей, Иван!")
# Пример 4: Преобразуване на низ в главни букви
upper_str <- toupper("здравей")
# Пример 5: Преобразуване на низ в малки букви
lower_str <- tolower("ЗДРАВЕЙ")

# Операция 10: Създаване и използване на фактори
# Пример 1: Създаване на фактор
factor_example <- factor(c("низък", "среден", "висок"))
# Пример 2: Достъп до нивата на фактор
levels(factor_example)
# Пример 3: Промяна на нивата на фактор
levels(factor_example) <- c("low", "medium", "high")
# Пример 4: Използване на фактор в анализ
summary(factor_example)

# Общо Задание: 
# Създайте вектор от числа от 1 до 10. Извършете следните операции:
# 1. Намерете средно аритметично на числата.
# 2. Създайте нов вектор, съдържащ само числата, които са по-големи от 5.
# 3. Преобразувайте втория вектор в фактор и изпишете неговите нива.
# 4. Създайте символен вектор с плодове и проверете дали "ябълка" е в него.

# Решение на Общото Задание:
# 1. Създаване на вектор
num_vector <- 1:10
# 2. Средно аритметично
mean_num_vector <- mean(num_vector)
# 3. Филтрация на вектора
filtered_vector <- num_vector[num_vector > 5]
# 4. Преобразуване в фактор
factor_vector <- factor(filtered_vector)
# 5. Показване на нивата на фактора
levels(factor_vector)
# 6. Създаване на символен вектор и проверка
fruit_vector <- c("ябълка", "банан", "портокал")
"ябълка" %in% fruit_vector


# Учебна Сесия 2: Импорт/Експорт на Данни и Основна Манипулация с Данни

# Операция 1: Инсталиране на пакети/библиотеки
# Пример 1: Стандартен подход - Инсталиране на ggplot2
install.packages("ggplot2")
# Зареждане на ggplot2 и прави достъпен за използване в текущата сесия на R
library(ggplot2)
# Пример 2: инсталиране на пакетите от "tidyverse", включително популярни пакети като ggplot2, dplyr, tidyr и други
install.packages("tidyverse")
# За да активиране на всички пакети от "tidyverse"
library(tidyverse)
# Пример 3: Инсталиране на специфичен пакет от github репозитор
# Инсталиране на devtools, ако не е инсталиран
install.packages("devtools")
# Зареждане на devtools
library(devtools)
# Инсталиране на библиотеката от GitHub "github_username/github_repository"
devtools::install_github("ramnathv/slidify")
# Инсталиране на ggplot2 от GitHub
devtools::install_github("tidyverse/ggplot2")

# Операция 2: Работа с данни от CSV файл
# Пример 1: Четене на CSV файл
df1 <- read.csv("mtcars.csv")
# Пример 2: Четене на CSV файл със заглавия
df2 <- read.csv("mtcars.csv", header = TRUE)
# Пример 3: Четене на CSV файл с разделител ';'
df3 <- read.csv("mtcars.csv", sep = ';')
# Пример 4: Запис на DataFrame в CSV
write.csv(df1, "output1.csv")
# Пример 5: Запис на DataFrame без редове с индекси
write.csv(df2, "output2.csv", row.names = FALSE)
# Пример 6: Запис на DataFrame с различен разделител
write.csv(df3, "output3.csv", sep = ';')

# Операция 3: Работа с данни от Excel файл
# Инсталиране на пакети readxl, openxlsx, ако все още не е инсталиран
install.packages("openxlsx")
install.packages("readxl")
# Зареждане на пакети readxl, openxlsx
library(openxlsx)
library(readxl)
# Пример 1: Четене на първия лист на Excel файл
df1 <- read_excel("01. Belgium.xlsx")
# Пример 2: Четене на конкретен лист със заглавия
# Заменете "SheetName" с името на листа
df2 <- read_excel("01. Belgium.xlsx", sheet = "Лист2")
df1 <- read_excel("01. Belgium.xlsx", sheet = "Лист1")
# Пример 3: Четене на Excel файл, като се задава броят на редовете за пропускане
# Пропуска първия ред
df3 <- read_excel("01. Belgium.xlsx", skip = 1)
# Пример 4: Запис на DataFrame в Excel файл
write.xlsx(df1, "output1.xlsx")
# Пример 5: Запис на DataFrame в различни листове на един Excel файл
write.xlsx(df1, "output1.xlsx", sheetName = "Sheet–1")
write.xlsx(df2, "output1.xlsx", sheetName = "Sheet–2")

# Операция 4: Работа със SPSS данни
# Пример 1: Четене на *.sav, *.zsav
# Инсталирайте и заредете пакета "haven"
install.packages("haven")
library(haven)
# прочетете SPSS файловете
data <- read_sav("mtcars.sav")
# Прегледайте първите няколко реда от данните
head(data)
# Пример 2: запис на данни във файл data.sav
write_sav(data, "data.sav")
# Пример за отваряне на SAS файл
sas_data <- read_sas("your_file.sas7bdat")

# Операция 5: Работа със .dbf данни
# Пример 1: Четене на *.dbf
# Инсталирайте и заредете пакета "foreign"
install.packages("foreign")
library(foreign)
# прочетете dbf файловете
data <- read.dbf("sample.dbf")
# Прегледайте първите няколко реда от данните
head(data)
# Пример 2: запис на DBF файл
write.dbf(data, "output.dbf")

# Операция 6: Използване на вградени набори от данни
# Пример 1: Използване на набора от данни mtcars
data(mtcars)
head(mtcars)
# Пример 2: Използване на набора от данни iris
data(iris)
head(iris)
# Пример 3: Използване на набора от данни airquality
data(airquality)
head(airquality)
# Пример 4: Извеждане на пълен списък от вградени данни
data()

# Операция 7: Създаване и използване на списъци
# Пример 1: Създаване на списък
list_example <- list(name = "Иван", age = 30, has_pet = TRUE)
# Пример 2: Достъп до елемент от списък
person_name <- list_example$name
# Пример 3: Промяна на елемент от списък
list_example$age <- 31

# Операция 8: индексация:
# Пример 1: Индексация на Вектор
vector_example <- c("a", "b", "c")
vector_example[2] # Връща "b"
vector_example[1:2] # Връща "a" и "b"
# Пример 2: Индексация на Матрица
matrix_example <- matrix(1:6, nrow = 2)
matrix_example[1, 2] # Връща елемента на първия ред и втората колона (2)
matrix_example[, 1] # Връща първата колона (1 3)
# Пример 3: Индексация на Лист
list_example <- list(letters, numbers)
list_example[[2]][3] # Връща третия елемент от втория лист (3)
# Пример 4: Индексация на Речник (list)
dictionary_example <- list(name = "John", age = 30)
dictionary_example$name # Връща стойността на ключ "name" (John)

# Операция 9: `[[ ]]` индексация за извлечане на стойността на елемента:
# Пример 1: Извличане на стойност от списък (list):
my_list <- list(a = 1, b = 2, c = 3)
value <- my_list[["b"]] # Връща стойността на ключ "b" (2)
# Пример 2: Извличане на стойност от колона в DataFrame:
df <- data.frame(Name = c("John", "Alice"), Age = c(30, 25))
ages <- df[["Age"]] # Връща стойностите на колоната "Age" (30, 25)
# Пример 3: Извличане на стойност от клетка в матрица (която може също да бъде част от DataFrame):
matrix_data <- matrix(1:6, nrow = 2)
value <- matrix_data[[2, 3]] # Връща стойността на втория ред и третата колона (6)

# Операция 10: Отпечатване на информация в конзолата
# Пример 1: Използване на print за отпечатване на информация с автоматично добавяне на нов ред
print("Това е съобщение, изпечатано с print.")
# Пример 2: Използване на cat за отпечатване на текст без автоматично добавяне на нов ред
cat("Това е", "съобщение,", "изпечатано с cat.")
cat("\n")  # Този ред добавя нов ред след отпечатването на текста
cat("Текстът след новия ред.")




# Общо Задание:
# 1. Създайте нов списък на име "person_info", който съдържа информация за човек, включваща следните ключове и стойности:
# - "first_name" - вашето име
# - "last_name" - вашата фамилия
# - "age" - вашата възраст
# - "city" - градът, в който живеете
# 2. Извличане и отпечатване с пояснения
# -име и фамилия на човека 
# -възрастта на човека
# -града, в който живее човекът

# Създаване на списъка "person_info"
person_info <- list(
  first_name = "Иван",
  last_name = "Петров",
  age = 30,
  city = "София"
)

# Извличане на името на човека
person_name <- person_info$first_name
cat("Име: ", person_info$first_name, person_info$last_name"\n","Възраст: ", person_age, "\n","Град: ", person_city, "\n")


# Учебна Сесия 3: Напреднала Манипулация с Данни (Data Frame)
data("mtcars")
df1 <- mcars

# Операция 1: създаване на DataFrame:
# Пример 1: Създаване на DataFrame от вектори
df1 <- data.frame(Name = c("John", "Alice"), Age = c(30, 25))
# Пример 2: Създаване на DataFrame от матрица
matrix_data <- matrix(1:6, nrow = 2)
df2 <- as.data.frame(matrix_data)
# Пример 3: Създаване на празен DataFrame с имена на колоните
df3 <- data.frame(Name = character(0), Age = numeric(0))

# Операция 2: извикване на елементи от DataFrame:
# Пример 1: Извикване на цяла колона
age_column <- df1$Age
# Пример 2: Извикване на отделен елемент по ред и колона
element <- df2[2, 3] # Връща елемента на втория ред и третата колона (6)
# Пример 3: Извикване на отделен елемент чрез $ оператор
name <- df1$Name[1] # Връща името на първия ред (John)

# Операция 7: Промяна на структурата на DataFrame
# Пример 1: Премахване на колона
df1$carb <- NULL
# Пример 2: Преименуване на колона
names(df1)[names(df1) == "X"] <- "model"
# Пример 3: Преорганизиране на колоните
df1 <- df1[c("model", "hp", "wt")]

# Операция 3: Създаване на нови колони в DataFrame
# Пример 1: Създаване на нова колона
df1$new_col <- df1$mpg * 2
# Пример 2: Създаване на колона с условие
df2$flag <- df2$mpg > 20
# Пример 3: Комбиниране на колони
df4 <- df2
df4$combined <- paste(df1$X, df2$hp)
# Пример 4: Създаване на деривативна колона
df1$ratio <- df1$hp / df1$drat

# Операция 4: Филтриране на редове в DataFrame
# Пример 1: Филтриране по числова стойност
filtered_df1 <- df1[df1$mpg > 20, ]
# Пример 2: Филтриране по символен низ
filtered_df2 <- df2[grepl("^M", df2$X), ]
# Пример 3: Комбинирано филтриране
filtered_df3 <- df1[(df1$mpg > 20) & grepl("^M", df2$X), ]

# Операция 5: Сортиране на данни в DataFrame
# Пример 1: Сортиране по една колона възходящо
sorted_df1 <- df1[order(df1$mpg), ]
# Пример 2: Сортиране по една колона низходящо
sorted_df2 <- df2[order(df2$hp, decreasing = TRUE), ]
# Пример 3: Сортиране по множество колони
sorted_df3 <- df1[order(df1$mpg, df1$hp), ]

# Операция 6: Агрегация на данни
# Пример 1: Изчисляване на средно аритметично
mean_value <- mean(df1$mpg)
# Пример 2: Намиране на максимална стойност
max_value <- max(df2$hp)
# Пример 3: Групиране и агрегация
aggreg_hp <-aggregate(hp ~ cyl, data = df1, mean)

# Операция 8: Съединяване на данни от различни DataFrames
# Пример 1: Съединяване по редове
combined_df1 <- rbind(df1, df1)
# Пример 2: Съединяване по колони
combined_df2 <- cbind(df1$X, df2$mpg)
# Пример 3: Съединяване с условие (JOIN)
combined_df3 <- merge(df1, df2, by = "X")
# Пример 4: Съединяване само на отделни колони - изисква library(dplyr)
combined_df4 <- merge(select(df1, X, mpg), select(df2, X, disp, hp), by = "X")

# Операция 12: Работа с фактори
# Пример 1: Преобразуване на числова колона във фактор
data(mtcars)
df1$factor_col <- factor(df1$cyl)
class(df1$factor_col)
class(df1$cyl)
# Пример 3: Преобразуване на символен низ във фактор
df2$flag <- df2$mpg > 20
df2$category <- as.factor(df2$flag)
# Пример 2: Промяна на нивата на фактор
levels(df1$factor_col) <- c("12", "10", "8")
# Пример 3: Сортиране на DataFrame според фактор
df1_sorted <- df1[order(df1$factor_col), ]

# Операция 10: Работа с дати и време
# Пример 1: Зареждане на набора от данни airquality
data("airquality")
df1 <- airquality
# Пример 2: Създаване на последователност от дати, съответстваща на дните в набора от данни airquality
# Предполага се, че данните започват от 1-ви май 1973 г.
dates <- seq(as.Date("1973-05-01"), by="day", length.out=nrow(airquality))
# Преобразуване на датите във формат на низ и добавяне към набора от данни airquality
df1$Dates <- as.Date(dates, format = "%Y-%m-%d")
# Пример 3: Извличане на годината от дата
df1$year <- format(df1$Dates, "%Y") 
# Пример 4: Изчисляване на разлика в дни между дати
days_diff <- as.numeric(difftime(df1$Dates[1], df1$Dates[length(df1$Dates)], units = "days"))
# Пример 5: Изчисляване на възраст в цели години - изисква library(lubridate)
age <- as.numeric(difftime(today(),as.Date("1953-11-28"))/365.2425)

df$date <- as.Date(paste(df$year, df$month, df$day, sep="-"), format="%Y-%m-%d")



# Операция 13: Работа с липсващи стойности
data("airquality")
df1 <- airquality
# Пример 1: Проверка за липсващи стойности
missing_values <- is.na(df1$Ozone)
# Пример 2: Замяна на липсващи стойности
df1$Ozone[is.na(df1$Ozone)] <- 0
df1$Ozone[is.na(df1$Ozone)] <-mean(df1$Ozone, na.rm = TRUE)
# Пример 3: Премахване на редове с липсващи стойности
df2 <- na.omit(df1)

# Общо Задание 1:
# Използвайте набора от данни mtcars. Извършете следните операции:
# 1. Създайте нова колона, която да представлява отношението на mpg към cyl.
# 2. Филтрирайте данните за автомобили с повече от 4 цилиндъра.
# 3. Сортирайте резултатите по новосъздадената колона в низходящ ред.
# 4. Изчислете средната стойност на hp за филтрираните автомобили.

# Решение на Общото Задание:
data(mtcars)
# 1. Създаване на нова колона
mtcars$ratio <- mtcars$mpg / mtcars$cyl
# 2. Филтриране на автомобили
filtered_cars <- mtcars[mtcars$cyl > 4, ]
# 3. Сортиране на данните
sorted_cars <- filtered_cars[order(filtered_cars$ratio, decreasing = TRUE), ]
# 4. Изчисляване на средната стойност на hp
mean_hp <- mean(sorted_cars$hp)

# Общо Задание 2:
# Използвайте набора от данни mtcars. Извършете следните операции:
# 1. Преобразувайте колоната 'mpg' във фактор, разделяйки стойностите на "Ниско", "Средно" и "Високо" спрямо тяхната стойност.
# 2. Изчислете средната мощност (hp) на автомобилите с "Високо" 'mpg'.
# 3. Създайте нов DataFrame, който съдържа само автомобили с "Ниско" и "Средно" 'mpg' и изчислете тяхната средна мощност.

# Решение на Общото Задание:
data(mtcars)
# 1. Преобразуване на 'mpg' във фактор
mtcars$mpg_factor <- cut(mtcars$mpg, breaks = c(-Inf, 15, 25, Inf), labels = c("Ниско", "Средно", "Високо"))
# 2. Средна мощност за 'Високо' mpg
high_mpg_hp <- mean(mtcars$hp[mtcars$mpg_factor == "Високо"])
# 3. Нов DataFrame и средна мощност
low_med_mpg <- mtcars[mtcars$mpg_factor %in% c("Ниско", "Средно"), ]
low_med_mpg_hp <- mean(low_med_mpg$hp)

                          
# Учебна Сесия 4: Основи на Автоматизацията в R                          

# Операция 1: Използване на функции от пакета dplyr
# `dplyr` е мощен и популярен пакет в R, който се използва за манипулация на данни. 
install.packages("dplyr") # освен ако не е вече инсталиран заедно с tidyvrse
library(dplyr)
data(mtcars)
# Пример 1: Филтриране с dplyr
filtered_df <- filter(df1, mpg > 20)
# Пример 2: Избиране на колони с dplyr
selected_df <- select(df2, cyl, 6)
# Пример 3: Сумаризиране на данни с dplyr
summary_df <- df1 %>% group_by(cyl) %>% summarise(Avg = mean(hp))
# Пример 4: Създаване на нови колони:
df_mutated <- mutate(mtcars, km_per_l = mpg * 0.425)
# Пример 5: Сортиране на редове:
df_sorted <- arrange(mtcars, desc(mpg))
# Пример 6: Групиране и агрегиране:
df_grouped <- group_by(mtcars, cyl)
summary_df <- summarise(df_grouped, avg_mpg = mean(mpg))
# Пример 7: Извличане на уникални редове:
df_unique <- distinct(mtcars, gear)
# Пример 8: Преименуване на колони:
df_renamed <- rename(mtcars, KilometersPerLiter = km_per_l)

# Операция 2: Работа с векторизирани операции
# Пример 1: Прилагане на функция по редове
df1$mean_col <- apply(df1[, c("hp", "mpg")], 1, sum)
# Пример 2: Прилагане на функция по колони
df1$mean_col <- apply(df1[, c("hp", "mpg")], 2, mean)
mean_col[1] <- mean(df1$Ozone, na.rm = TRUE)
# Пример 1: Използване на lapply
df1$new_col <- lapply(df1$hp, function(x) x^2)
# Пример 2: Използване на sapply
df1$new_col <- sapply(df1$mpg, sqrt)
# Пример 3: Използване на mapply
df1$new_col <- mapply(function(x, y) x * y, df1$cyl, df1$hp)
                      
# Операция 4: Използване на математически функции
# Пример 1: Корен квадратен
sqrt_vector <- sqrt(vector1)
# Пример 2: Логаритъм
log_vector <- log(vector1)
# Пример 3: Експоненциална функция
exp_vector <- exp(vector1)

# Операция 5: Използване на условни оператори
# Пример 1: if конструкция
if (addition > 7) {
  result <- "По-голямо от 7"
}
# Пример 2: if-else конструкция
if (subtraction > 0) {
  result <- "Положително"
} else {
  result <- "Отрицателно или нула"
}
# Пример 3: if-else if-else конструкция
if (multiplication < 10) {
  result <- "По-малко от 10"
} else if (multiplication == 15) {
  result <- "Равно на 15"
} else {
  result <- "По-голямо от 15"
}
# Пример 4: ifelse за създаване на нова колона
data("mtcars")
df1 <- mtcars
df1$new_col <- ifelse(df1$hp > 100, "Sports", "Econnomic")

# Операция 6: Използване на цикли
# Пример 1: for цикъл
for (i in 1:5) {
  print(i)
}
# Пример 2: for цикъл за обработка на данни
for (i in 1:nrow(df1)) {
  df1$new_col[i] <- df1$cyl[i] * df1$hp[i]
}
# Пример 3: while цикъл
i <- 1
while (i <= 5) {
  print(i)
  i <- i + 1
}
# Пример 4: while цикъл за обработка на данни
i <- 1
while (i <= 10) {
  df1$new_col[i] <- df1$cyl[i] / df1$hp[i]
  i <- i + 1
}
# Пример 5: repeat цикъл с условие за прекъсване
i <- 1
repeat {
  print(i)
  if (i >= 5) {
    break
  }
  i <- i + 1
}

# Операция 7: Създаване на потребителски функции
# Пример 1: Функция без аргументи
say_hello <- function() {
  print("Здравей!")
}
# Пример 2: Функция с един аргумент
square <- function(x) {
  x^2
}
# Пример 3: Функция с множество аргументи
divide <- function(x, y=10) {
  x / y
}


# Обшо Задание 1: Използвайте пакета `dplyr`, за да изпълните следните операции върху данните от `iris`:
# Филтрирайте редовете на `iris`, така че да включват само наблюдения със средна дължина на чашката (`Sepal.Length`) по-голяма от 5 и средна ширина на чашката (`Sepal.Width`) по-малка от 3.
# Създайте нова колона с име "Petal.Area", която представлява произведението на дължината на листото (`Petal.Length`) и ширината на листото (`Petal.Width`) за всяко наблюдение.
# Групирайте данните по вид цвете (`Species`) и намерете средната стойност на "Petal.Area" за всяка група.
# Сортирайте групите по средната стойност на "Petal.Area" в намаляващ ред.
# Изберете само първите 5 групи (с най-голяма средна стойност на "Petal.Area") и запазете резултата в нова променлива.

# Решение
# Използване на пакета dplyr
library(dplyr)
# Филтриране на данните
filtered_iris <- iris %>%
filter(Sepal.Length > 5, Sepal.Width < 3)
# Създаване на нова колона "Petal.Area"
mutated_iris <- filtered_iris %>%
mutate(Petal.Area = Petal.Length * Petal.Width)
# Групиране по вид цвете и изчисляване на средната стойност на "Petal.Area" за всяка група
grouped_iris <- mutated_iris %>%
group_by(Species) %>%
summarise(Avg_Petal.Area = mean(Petal.Area))
# Сортиране на групите в намаляващ ред по средната стойност на "Petal.Area"
sorted_iris <- grouped_iris %>%
arrange(desc(Avg_Petal.Area))
# Избор на първите 5 групи
top_5_species <- head(sorted_iris, 5)

# Общо Задание 2:
# Напишете функция, която приема числов вектор и връща нов вектор, където всеки елемент е квадратът на оригиналния, ако е по-голям от 10, и оригиналния елемент в противен случай.
# Тествайте функцията с вектор от числата 1 до 20.

# Решение на Общото Задание:
square_if_gt_10 <- function(vector) {
  sapply(vector, function(x) ifelse(x > 10, x^2, x))
}
test_vector <- 1:20
result_vector <- square_if_gt_10(test_vector)
print(result_vector)

         
# Учебна Сесия 5: Въведение в Визуализацията на Данни

# Операция 1: Създаване на основни графики
data("mtcars")
# Пример 1: Хистограма
hist(mtcars$mpg)
# Пример 2: Точков диаграма
plot(mtcars$mpg, mtcars$hp)
# Пример 3: Линейна диаграма
plot(mtcars$mpg, type = "l")

# Операция 2: Персонализиране на графики
# Пример 1: Добавяне на заглавия и означения на оста
plot(mtcars$mpg, main = "Разпределение на MPG", xlab = "MPG", ylab = "Честота")
# Пример 2: Промяна на цвета и типа на точките
plot(mtcars$mpg, mtcars$hp, col = "red", pch = 19)
# Пример 3: Добавяне на линии на решетката
plot(mtcars$mpg, mtcars$hp); grid()

# Операция 3: Създаване на сложни графики
# Пример 1: Стекова стълбова диаграма
barplot(table(mtcars$cyl, mtcars$gear), col = c("red", "blue", "green"), legend = TRUE)
# Пример 2: Кутия за уси (Boxplot)
boxplot(mtcars$mpg ~ mtcars$cyl)
# Пример 3: Графика на разсеяване със фасети
ggplot(mtcars, aes(x = mpg, y = hp)) + geom_point() + facet_wrap(~cyl)

# Операция 4: Визуализация на времеви редове
# Пример 1: Линейна диаграма на времеви ред
plot(AirPassengers, type = "l")
# Пример 2: Декомпозиция на времеви ред
library(stats)
plot(decompose(AirPassengers))
# Пример 3: Автокорелация на времеви ред
acf(AirPassengers)

# Операция 5: Интерактивни графики
# Пример 1: Интерактивна точкова диаграма
# install.packages("plotly")
library(plotly)
ggplotly(ggplot(mtcars, aes(x = mpg, y = hp)) + geom_point())

# Операция 7: Визуализация на категорийни данни
# Пример 1: Стълбова диаграма
barplot(table(mtcars$cyl))
# Пример 2: Торта диаграма
pie(table(mtcars$cyl))
# Пример 3: Хийтмап
heatmap(as.matrix(mtcars[, 1:4]))

# Операция 6: Запазване на графики като изображения
data(mtcars)
# Пример 1: Запазване на графика като PNG
png("my_plot.png")
plot(mtcars$mpg, mtcars$hp)
dev.off()
# Пример 2: Запазване на графика като PDF
pdf("my_plot.pdf")
plot(mtcars$mpg, mtcars$hp)
dev.off()
# Пример 3: Запазване на графика като JPEG
jpeg("my_plot.jpeg")
plot(mtcars$mpg, mtcars$hp)
dev.off()
# Пример 4: Запис в SVG файл
svg("myplot.svg")
plot(mtcars$mpg)
dev.off()
         
# Операция 8: Работа с knitr за създаване на доклади
install.packages("knitr")
library(knitr)
# Създайте файл с разширение `.Rmd` и добавете следното съдържание:
---
title: "Примерен Доклад"
author: "Вашето Име"
date: "Дата"
output: html_document
---

Това е примерен доклад.

#За да компилирате R Markdown документ - бутона "Knit" в RStudio или
rmarkdown::render("your_document.Rmd")
         
# Пример 1: Вграждане на R код в Markdown документ
```{r}
summary(cars)
# Пример 2: Графика:
```{r pressure, echo=FALSE}
plot(pressure)
# Пример 3: Добавяне на описания към кода в Markdown
```{r, echo=TRUE}
# # Това е графика
plot(pressure)
# Пример 4: Контролиране на размера на изображението в Markdown
```{r, fig.height=8, fig.width=10}
plot(pressure)
# Пример 5: Добавяне на Опции към R Кодовите Блокове
# Можете да контролирате поведението на всеки R кодов блок с опции:
{r results='hide'}
# Кодът се изпълнява, но резултатите не се показват в крайния документ
summary(cars)
{r echo=FALSE}
# Кодът се изпълнява и резултатите се показват, но самият код не се показва
plot(cars)
# Пример 6: Включване на LaTeX за Създаване на PDF Документ
#Ако искате да създадете PDF, можете да включите LaTeX команди в R Markdown документа:

---
title: "Примерен Доклад"
output: 
  pdf_document:
    latex_engine: xelatex
---
  
  Това е примерен доклад.

Можете да пишете LaTeX формули като:
  
  $$ e^{i\pi} + 1 = 0 $$
  
  ```{r}
# Примерен R код
plot(cars)

         
# Общо Задание:
# Използвайте набора от данни mtcars. Извършете следните операции:
# 1. Създайте хистограма на 'mpg'.
# 2. Създайте точкова диаграма на 'mpg' спрямо 'hp', обозначете различните цилиндри ('cyl') с различни цветове.
# 3. Създайте линейна диаграма на 'mpg' спрямо 'wt', добавете линия на тренда.
# 4. Публикувайте всички графики в общ доклад, използвайки knittr.

# Решение на Общото Задание:
---
title: "Доклад за mtcars"
output: html_document
---

# Графики и Анализ

```{r}
# Графика 1: Хистограма на 'mpg'
hist(mtcars$mpg, main = "Разпределение на MPG", xlab = "MPG")
# Графика 2: Точкова диаграма на 'mpg' спрямо 'hp'
plot(mtcars$mpg, mtcars$hp, col = mtcars$cyl, pch = 19, main = "MPG спрямо HP", xlab = "MPG", ylab = "HP")
legend("topright", legend = unique(mtcars$cyl), col = unique(mtcars$cyl), pch = 19)
# Графика 3: Линейна диаграма на 'mpg' спрямо 'wt' с линия на тренда
plot(mtcars$mpg, mtcars$wt, main = "MPG спрямо WT", xlab = "MPG", ylab = "WT")
abline(lm(mtcars$wt ~ mtcars$mpg), col = "red")


# Учебна Сесия 6: Напреднала Визуализация с ggplot2 и други

# Операция 1: Работа с ggplot2
install.packages("ggplot2")
library(ggplot2)
#1. Основна Структура:
#`ggplot2` използва граматика за графики, която включва определение на данни, естетика (аesthetic mappings) и слоеве.
# ggplot(data = <данни>, aes(<естетика>)) + <слой>
#2. Данни и Естетика (aes):
# data =
# aes(x = <x-променлива>, y = <y-променлива>)
#3. Слоевете се добавят с `+`.
# geom_point() + geom_line()
#4. Скали и означения (labels)и др.
# scale_x_continuous(name = "X axis")
#5. Разделяне на данните в множество панели.
# facet_wrap(~ <фактор>) или facet_grid(<рeдове> ~ <колони>)
#6. Визуален стил
# theme_minimal(), theme_classic() и т.н.
#7. Заглавия и Етикети.
# labs(title = "<заглавие>", x = "<означение x-ост>", y = "<означение y-ост>")
#8. Логаритмични трансформации на остите.
# scale_y_log10()
#9. Превърнете `ggplot` графика в интерактивна с `ggplotly()` от пакета `plotly`.
#10. Запазете графиката във файл с `
# ggsave("filename.png")
# Пример 1: Основна Точкова Диаграма (Scatter Plot)
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point()
# Пример 2: Линейна Диаграма (Line Plot)
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_line() +
  geom_point()
# Пример 3: Персонализиране с теми
ggplot(mtcars, aes(x = mpg, y = hp)) + geom_point() + theme_minimal()

# Операция 2: Други диаграми с ggplot2
# Пример 1: Хистограма (Histogram)
ggplot(mtcars, aes(x = mpg)) +
  geom_histogram(bins = 15, fill = "blue", color = "black")
# Пример 2: Кутиева Диаграма (Boxplot)
ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_boxplot()
# Пример 3: Bar Plot с Групиране
ggplot(mtcars, aes(x = factor(cyl), fill = factor(gear))) +
  geom_bar(position = "dodge")
# Пример 4: Фасетна Диаграма (Facet Plot)
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  facet_wrap(~cyl)

# Операция 3: Смесени диаграми с ggplot2
# Пример 1: Графика с Гладки Линии (Smooth Line Plot)
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_smooth(method = "loess", se = FALSE)
# Пример 2: Създаване на графика с множество графични елементи
ggplot(mtcars, aes(x = mpg, y = hp, color = factor(cyl))) + geom_point() + geom_line()
# Пример 3: Графика с Наложени Хистограми (Overlaid Histograms)
ggplot(mtcars, aes(x = mpg, fill = factor(cyl))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 15)
# Пример 4: Графика с Етикети (Plot with Labels)
ggplot(mtcars, aes(x = wt, y = mpg, label = rownames(mtcars))) +
  geom_text(check_overlap = TRUE, vjust = 1.5)
# Пример 5: Интерактивна Графика (Interactive Plot) с plotly
library(plotly)
p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
ggplotly(p)

# Операция 4: Използване на различни цветови палитри
# Пример 1: Използване на вградени цветови палитри
ggplot(mtcars, aes(x = mpg, y = hp, color = factor(cyl))) + geom_point() + scale_color_brewer(palette = "Dark2")
# Пример 2: Дефиниране на собствена цветова палитра
my_colors <- c("red", "green", "blue")
ggplot(mtcars, aes(x = mpg, y = hp, color = factor(cyl))) + geom_point() + scale_color_manual(values = my_colors)
# Пример 3: Използване на палитра от RColorBrewer
library(RColorBrewer)
display.brewer.all()

# Операция 5: Визуализация на комплексни структури от данни
# Пример 1. Топлинна Карта (Heatmap)
ggplot(mtcars, aes(x = factor(cyl), y = factor(gear))) +
  geom_tile(aes(fill = mpg))
# Пример 2. Графика с Дървена Структура (Dendrogram)
install.packages("ggdendro")
library(ggdendro)
hc <- hclust(dist(mtcars))
dendro <- as.dendrogram(hc)
ggdendrogram(dendro)
# Пример 3. Корелационна Матрица (Correlation Matrix)
install.packages("reshape2")
library(reshape2)
corr_matrix <- cor(mtcars)
melted_corr <- melt(corr_matrix)
ggplot(melted_corr, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile()
# Пример 4: 3D точкова графика
install.packages("scatterplot3d")
library(scatterplot3d)
scatterplot3d(mtcars$mpg, mtcars$hp, mtcars$cyl)
         
# Операция 5: Графика с Географски Данни (Geographical Plot)
install.packages("maps")
library(maps)
# Създаване на набор от данни за Европа със стойности за всяка държава
europe_data <- map_data("world") %>%
  filter(tolower(region) %in% c("france", "spain", "germany", "italy", "poland", "norway", "sweden", "finland", "ukraine", "greece", "bulgaria", "romania"))
europe_data$values <- runif(n = nrow(europe_data), min = 1, max = 100)
# Създаване на картата
ggplot(europe_data, aes(x = long, y = lat, group = group, fill = values)) +
  geom_polygon(color = "black") +
  theme_minimal() +
  scale_fill_viridis_c() + # използване на цветова скала
  labs(fill = "Стойност", title = "Карта на Европа") +
  coord_fixed(1.3) # фиксиране на съотношението на аспектите

# Пример 6. Мрежа от Графики (Network Graph)
install.packages("igraph")
library(igraph)
# Създаване на основния пръстен
g <- make_ring(10)
# Добавяне на допълнителни върхове
g <- add_vertices(g, nv = 10, color = c(rep("blue", 3), rep("green", 7)))
# Добавяне на допълнителни ребра
g <- add_edges(g, c(1, 11, 11, 12, 12, 13, 13, 1))
g <- add_edges(g, c(2, 14, 14, 15, 15, 16, 16, 2))
g <- add_edges(g, c(3, 17, 17, 18, 18, 3))
# Определяне на размерите на върховете
V(g)$size <- rep(5, vcount(g))  # Задаване на размер 5 за всички върхове
V(g)$size[11:20] <- 10          # Промяна на размера на новодобавените върхове на 10
# Определяне на цвета на върховете
V(g)$color <- rep("gray", vcount(g))
V(g)$color[11:20] <- "red"
# Определяне на цвета на ребрата
E(g)$color <- "gray"
# Визуализация на графиката
plot(g, vertex.label = NA, edge.arrow.size = 0.5, vertex.size = V(g)$size,
     vertex.color = V(g)$color, edge.color = E(g)$color)

# Пример 7. Визуализация на Текстови Данни (Word Cloud)
install.packages("wordcloud")
library(wordcloud)
install.packages("tm")
library(tm)
install.packages("rvest")
library(rvest)
# Изтегляне на HTML съдържанието на статията
text <- tolower(readLines("text.txt", warn = FALSE))
text <- paste(text, collapse = " ")
# Създаване на корпус от текст
corpus <- Corpus(VectorSource(text))
corpus <- tm_map(corpus, removeWords, stopwords("english"))
dtm <- DocumentTermMatrix(corpus)
# Превръщане на матрицата в обект с думи и техните честоти
word_freq <- as.data.frame(as.table(dtm))
# Създаване на облак от думи
# Преобразуване на колоната "Frequency" към числа
word_freq$Frequency <- as.numeric(word_freq$Freq)
# Създаване на облак от думи
wordcloud(words = word_freq$Terms, freq = word_freq$Frequency)

# Пример 8. Радарова Графика (Radar Chart)
install.packages("fmsb")
library(fmsb)
data("mtcars")
# Подготовка на данните за радарната диаграма
mtcars_scaled <- as.data.frame(lapply(mtcars[1:5,], function(x) {
  (x - min(x)) / (max(x) - min(x))
}))
# Preserve the original row names (car models) after scaling
rownames(mtcars_scaled) <- rownames(mtcars[1:5, ])
# Create a row for min and max values
min_row <- setNames(rep(0, ncol(mtcars_scaled)), names(mtcars_scaled))
max_row <- setNames(rep(1, ncol(mtcars_scaled)), names(mtcars_scaled))
mtcars_scaled <- rbind(max_row, min_row, mtcars_scaled)
rownames(mtcars_scaled)[1] <- "Max"
rownames(mtcars_scaled)[2] <- "Min"
# Създаване на радарна диаграма
radarchart(mtcars_scaled)

# Общо Задание:
# Използвайте набора от данни mtcars. Създайте интерактивна графика на 'mpg' спрямо 'hp' с plotly, където различните цилиндри ('cyl') са обозначени с различни цветове. Запазете графиката като HTML файл.

# Решение на Общото Задание:
library(plotly)
p <- ggplot(mtcars, aes(x = mpg, y = hp, color = factor(cyl))) + geom_point()
interactive_plot <- ggplotly(p)
htmlwidgets::saveWidget(interactive_plot, "mtcars_plot.html")

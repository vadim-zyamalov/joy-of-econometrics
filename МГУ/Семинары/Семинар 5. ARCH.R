# Библиотека Quantmod дает нам удобый инструмент для загрузки данных из различных источников, но только западных.
library(magrittr)
library(quantmod)
library(urca)
library(dynlm)

# Загрузим дневные данные по торгам акциями "одной фруктовой компании".
getSymbols("AAPL", from = "2020/01/01", to = "2025/12/31", periodicity = "daily")

# Посчитаем логарифм цены
apple_l <- log(AAPL$AAPL.Adjusted)
lp <- zoo(apple_l)

# Посчитаем логарифм доходности...
apple_r <- log(AAPL$AAPL.Adjusted / Lag(AAPL$AAPL.Adjusted, 1))
# и отсечем первое -- NA -- значение.
apple_r <- na.omit(apple_r)
# names(apple_r) <- "lr"

N <- length(lp)

# Процедура Доладо-Дженкинсона-Сосвилья-Риверо
ur.df(apple_l, type = "trend", lags = N^(1 / 3), selectlags = "BIC") |> summary()

dynlm(d(lp) ~ trend(lp) + L(d(lp), 1:1)) |> summary()

ur.df(apple_l, type = "drift", lags = N^(1 / 3), selectlags = "BIC") |> summary()

dynlm(d(lp) ~ L(d(lp), 1:1)) |> summary()

ur.df(apple_l, type = "none", lags = N^(1 / 3), selectlags = "BIC") |> summary()

# Посмотрим на график.
plot.ts(apple_r)

# Удостоверимся, что ряд доходностей стационарен.
ur.df(apple_r, type = "drift", lags = 10, selectlags = "AIC") |> summary()
ur.pp(apple_r, type = "Z-tau", model = "constant", lags = "long") |> summary()
ur.kpss(apple_r, type = "mu", lags = "long") |> summary()

# Построим графики ACF и PACF.
acf(apple_r)
pacf(apple_r)

# При помощи функции auto.arima можно автоматически подобрать "оптимальную" SARIMAX модель.
# В нашем случае все ограничится ARMA :)
model1 <- forecast::auto.arima(apple_r)
summary(model1)

# Построим ряд квадратов остатков...
residsq <- model1$residuals^2
# и посмотрим на графики ACF и PACF
acf(residsq)
pacf(residsq)

# Если остатки получены не из ARMA модели, то надо превратить их во временной ряд.
# residsq <- zoo(residsq)

# Оценим AR(4) для квадратов остатков.
# Лаги создаются не по-одиночке, а разом, при помощи последовательности 1:4.
model2 <- dynlm(residsq ~ L(residsq, 1:10))
summary(model2)

# Протестируем совместную гипотезу о равенстве 0 коэффициетов при лагах.
# Гипотезы генерируются при помощи функции paste.
# Вообще, можете вписывать гипотезы вручную в виде:
# с("??? = ?", "??? = ?", ...)
car::linearHypothesis(model2, paste("L(residsq, 1:10)", 1:10, "=0", sep = ""))
# Либо проведите тест Вальда: первый параметр --- ковариационная матрица, второй --- оценки коэффициентов,
# третий --- номера зануляемых коэффициентов.
aod::wald.test(Sigma = vcov(model2), b = coef(model2), Terms = 2:7)

# Протестируем первую модель на ARCH при помощи специальной функции.
aTSA::arch.test(model1)

.model1 <- model1
class(.model1) <- "Arima"
aTSA::arch.test(.model1)

# Оценим GARCH(1,1) при помощи пакета `rugarch`.
# Сначала надо создать объект с описанием модели. В следующем примере создется модель типа
# ARMA(0,1)-GARCH(1,1).
rugarch::ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 1)
  )
) |> rugarch::ugarchfit(data = apple_r)

# Следующая модель имеет спецификацию
# ARMA(0,0)-GARCH(1,1).
rugarch::ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0)
  )
) |> rugarch::ugarchfit(data = apple_r)

# Тест на нормальность
fBasics::lillieTest(apple_r)

# Следующая модель имеет спецификацию
# ARMA(0,0)-GARCH(1,1).
# Используется распределение Стьюдента.
rugarch::ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0)
  ),
  distribution.model = "std"
) |> rugarch::ugarchfit(data = apple_r)

# ARCH-M
rugarch::ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    archm = TRUE,
    archpow = 2
  ),
  distribution.model = "std"
) |> rugarch::ugarchfit(data = apple_r)

# IGARCH
rugarch::ugarchspec(
  variance.model = list(
    model = "iGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0)
  ),
  distribution.model = "std"
) |> rugarch::ugarchfit(data = apple_r)

# eGARCH
rugarch::ugarchspec(
  variance.model = list(
    model = "eGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0)
  ),
  distribution.model = "std"
) |> rugarch::ugarchfit(data = apple_r)

# Скачаем данные по акциям IBM.
getSymbols("IBM", from = "2020/01/01", to = "2025/12/31", periodicity = "daily")

ibm_r <- log(IBM$IBM.Adjusted / Lag(IBM$IBM.Adjusted, 1))
ibm_r <- na.omit(ibm_r)

plot.ts(ibm_r)

# Объединим доходности в единый массив.
# Так как оба ряда --- временные ряды (не потому, что это так, а потому, что они --- объекты xts!),
# их можно объединить при помощи функции merge.
# Затем заполним пропуски 0.
data_r <- merge(apple_r, ibm_r) |> na.fill(0)

# Эта функция работает долго (минут 10) и не сходится...
# Иронично, что пакет MTS написан для учебника по этим моделям...
# res_BEKK <- MTS::BEKK11(data_r)

# Следующий пакет создан для оценивания BEKK модели.
# Работает ОЧЕНЬ шустро.
# Сначала создается объект с описанием модели.
BEKKs::bekk_spec() |>
  # Затем модель оцениванится.
  BEKKs::bekk_fit(data_r) |> summary()

# В лекции 1 приведена модель CCC. Далее оценивается ее расширение.
# Сначала проводится предварительный этап --- оценивание отдельных моделей для рядов.
pre_DCC <- MTS::dccPre(data_r, cond.dist = "std")
# Потом на базе этих оценок оцениваются параметры полинома, описывающего условную корреляцию между рядами.
res_DCC <- MTS::dccFit(pre_DCC$sresi, cond.dist = "std")

# Альтернативная функция.
# Сначала из отдельных описаний уравнений (см. выше пакет rugarch) создается сборный объект.
mres_DCC <- rugarch::ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0)
  )
) %>%
  replicate(2, .) |>
  rugarch::multispec() |>
  # Затем на его основе строится спецификация DCC модкли.
  rmgarch::dccspec(distribution = "mvt", model = "DCC") |>
  # И, наконец, оценка!
  rmgarch::dccfit(data = data_r)
mres_DCC

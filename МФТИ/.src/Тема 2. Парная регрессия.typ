#import "@preview/theorion:0.6.0": *
#import "@preview/tschich:0.2.0": *
#import "@local/neat-document:0.1.0": *
#import "@local/trinkets:0.1.0": *
#set document(title: [Парная линейная регрессия])
#set table(stroke: none, inset: 0.5em)
#show: neat-document.with(margin: tschich-var(210mm, 297mm, 1 / 10))
#show: show-theorion

После тестирования гипотез о равенстве средних нас может заинтересовать вопрос: как изменится одна из переменных при изменении какой-либо другой?
Очевидным будет получение оценки следующего параметра
$ beta_1 = frac(Delta Y, Delta X) $
Нетрудно увидеть, что данный параметр является частью уравнения прямой
$ Y = beta_0 + beta_1 X $
где $beta_0 = Y_0 - beta_1 X_0$.

Этими несложными рассуждениями мы практически ввели то, что мы будем называть моделью *парной линейной регрессии*.
Для окончательного определения модели нам надо посмотреть на то, как данные выглядят в реальной жизни.

#align(center)[#image(height: 40%, "fig/likbez_3_1.png")]

На графике представлены сгенерированные, но вполне возможные в реальности данные и линия, "наилучшим" способом их описывающая.
Как можно увидеть, реальные данные не лежат идеально на построенной линии.
Тем самым можно предположить, что в данных, помимо связи между переменными, есть еще что-то.
Так как отклонения не несут систематического характера, то мы можем рассматривать их как *случайное возмущение*.

Таким образом, для каждого конкретного наблюдения уравнение прямой примет вид
#math.equation(block: true, numbering: "(1)")[$Y_i = beta_0 + beta_1 X_i + u_i$] <eq1>
где
/ $Y_i$: объясняемая переменная. Иногда, очень редко, ее называют _регрессанд_.
/ $X_i$: объясняющая переменная. Очень часто ее называют _регрессор_.
/ $u_i$: случайное возмущение, называемое обычно _случайной ошибкой_.

Само уравнение называется *моделью парной линейной регрессии* или *моделью наблюдений*.

Приведенное выше уравнение
$ Y = beta_0 + beta_1 X $
называется *линией теоретической регрессии* или *моделью связи*.

= Оценивание коэффициентов модели

Как же нам оценить параметры модели?

Рассмотрим модель связи.
В *идеале* все точки из нашей выборки должны удовлетворять этой модели, то есть должно выполняться
$ Y_i = beta_0 + beta_1 X_i, quad forall i $
Тем самым мы получаем систему $n$ уравнений c *двумя* неизвестными.
В общем случае такая система не имеет точного решения, однако возможно найти решение _приблизительное_.

Для этого посмотрим на невязку, равную $Y_i - beta_0 - beta_1 X_i, quad forall i$.
Попробуем подобрать значения параметров $beta_0$ и $beta_1$ так, чтобы некое агрегированное значение этих невязок было минимальным.
На практике часто используют сумму квадратов, получая следующую задачу
$
  hat(beta)_0, hat(beta)_1 = arg min_(hat(beta)_0, hat(beta)_1) sum_(i = 1)^n (Y_i - hat(beta)_0 - hat(beta)_1 X_i)^2 eq.colon arg min_(hat(beta)_0, hat(beta)_1) Q
$
Данная задача решается элементарно: достаточно взять частные производные по неизвестным параметрам, приравнять их нулю и решить полученную систему
$
  frac(partial Q, partial hat(beta)_0) & = - 2 sum_(i = 1)^n (Y_i - hat(beta)_0 - hat(beta)_1 X_i)     & = 0 \
  frac(partial Q, partial hat(beta)_1) & = - 2 sum_(i = 1)^n X_i (Y_i - hat(beta)_0 - hat(beta)_1 X_i) & = 0
$
Нетрудно получить, что
$
  hat(beta)_0 & = macron(Y) - hat(beta)_1 macron(X)\
  hat(beta)_1 & = frac(sum_(i = 1)^n (Y_i - macron(Y))(X_i - macron(X)), sum_(i = 1)^n (X_i - macron(X))^2) = s_(Y X) / s_X^2 = beta_1 + frac(sum_(i = 1)^n (X_i - macron(X))u_i, sum_(i = 1)^n (X_i - macron(X))^2)
$

После оценивания модели получим *линию выборочной регрессии*
$ Y_i = hat(beta)_0 + hat(beta)_1 X_i + hat(u)_i $
где $hat(u)_i$ --- *остаток* регрессии.

== Предпосылки МНК

Полученные выше оценки являются <<хорошими>> при выполнении следующих предпосылок.

0. Мы *действительно* имеем дело с линейной моделью!
1. $Ex (u_i vert X_i) = 0$.
  Данная предпосылка говорит о том, что наблюдения могут отклоняться от линии регрессии в моменте, но в среднем эти отклонения должны быть нулевыми.
  #align(center)[#image(height: 15%, "fig/OLS_homo_errors.pdf")]
2. $(Y_i, X_i) tilde.op "iid"$.
  Данная предпосылка означает, что выборка должны быть случайной и независимой.
3. $Y_i$ и $X_i$ должны иметь конечные четвертые моменты.
  Иначе это формулируется как отсутствие сильных выбросов.
  #align(center)[#image(height: 25%, "fig/OLS_outliers.pdf")]

При выполнении этих предпосылок,
+ $sqrt(n)(hat(beta)_1 - beta_1) tilde.op_"asy" N(0, sigma_(hat(beta)_1)^2)$ где
  $ sigma^2_(hat(beta)_1) = 1/n (Var [(X_i - Ex(X)] u_i]) / ([Var(X_i)]^2) $
+ $sqrt(n) (hat(beta)_0 - beta_0) tilde.op_"asy" N(0, sigma_(hat(beta)_0)^2)$ где
  $
    sigma^2_(hat(beta)_0) & = 1/n (Var(H_i u_i))/((Ex[H^2_i])^2) \
                      H_i & = 1 - (Ex(X))/(Ex(X_i^2)) X_i
  $

== Свойства МНК-оценок

+ $sum_(i = 1)^n hat(u)_i = 0$
  Данное свойство следует из условия первого порядка для $frac(partial Q, partial hat(beta)_0)$.
  Иначе его справедливость можно показать следующим образом:
  $
    sum_(i = 1)^n hat(u)_i & = sum_(i = 1)^n (Y_i - hat(beta)_0 - hat(beta)_1 X_i) = \
                           & = n macron(Y) - n hat(beta)_0 - n hat(beta)_1 macron(X) = \
                           & = n [macron(Y) - (macron(Y) - hat(beta)_1 macron(X)) - hat(beta)_1 macron(X)] = 0
  $

+ $sum_(i = 1)^n hat(u)_i X_i = 0$.
  Отсюда, а также из первого свойства, следует, что $s_(hat(u)X) = 0$.
  Само это свойство является вторым условием первого порядка для $frac(partial Q, partial hat(beta)_1)$.

+ $1 / n sum_(i = 1)^n hat(Y) = macron(Y)$.
  $
    1 / n sum_(i = 1)^n hat(Y) & = 1 / n sum_(i = 1)^n (hat(beta)_0 + hat(beta)_1 X_i) = \
                               & = hat(beta)_0 + hat(beta)_1 macron(X) = \
                               & = macron(Y) - hat(beta)_1 macron(X) + hat(beta)_1 macron(X) = macron(Y)
  $

+ $(macron(Y), macron(X)) in hat(Y) = hat(beta)_0 + hat(beta)_1 X$.

+ $italic("TSS") =italic("ESS") + italic("RSS")$, где
  - $italic("TSS") = sum_(i = 1)^n (Y_i - macron(Y))^2$ --- полная сумма квадратов (total sum of squares);
  - $italic("ESS") = sum_(i = 1)^n (hat(Y)_i - macron(Y))^2$ --- объясненная сумма квадратов (explained sum of squares);
  - $italic("RSS") = sum_(i = 1)^n (Y_i - hat(Y)_i)^2 = sum_(i = 1)^n hat(u)_i^2$ --- остаточная сумма квадратов или кумма квадратов остатков (residual sum of squares).
  $
    italic("TSS") = sum_(i = 1)^n (Y_i - macron(Y))^2 & = sum_(i = 1)^n [(Y_i - hat(Y)) + (hat(Y) - macron(Y))]^2 = \
    & = underbrace(sum_(i = 1)^n (Y_i - hat(Y))^2, italic("RSS")) + underbrace(sum_(i = 1)^n (hat(Y) - macron(Y))^2, italic("ESS")) + 2 sum_(i = 1)^n (Y_i - hat(Y))(hat(Y) - macron(Y)) = bigast \
    hat(Y)_i - macron(Y) & = hat(beta)_0 + hat(beta)_1 X_i - macron(Y) = \
    & = macron(Y) - hat(beta)_1 macron(X) + hat(beta)_1 X_i - macron(Y) = \
    & = hat(beta)_1 (X_i - macron(X)) \
    bigast & = italic("RSS") + italic("ESS") + 2 hat(beta)_1 sum_(i = 1)^n hat(u)_i (X_i - macron(X)) = \
    & = italic("RSS") + italic("ESS") + 2 hat(beta)_1 autolr(underbrace(sum_(i = 1)^n hat(u)_i X_i, = 0) - macron(X) underbrace(sum_(i = 1)^n hat(u)_i, = 0)) =\
    & = italic("RSS") + italic("ESS")
  $

== Графическая интерпретация

Представим @eq1 в следующем виде
$
  underbrace(vec(y_1, dots.v, y_N), Vect(y)) = beta_0 underbrace(vec(1, dots.v, 1), Vect(1)) + beta_1 underbrace(vec(x_1, dots.v, x_N), Vect(x)) + underbrace(vec(u_1, dots.v, u_N), Vect(u))
$
то есть в виде линейной комбинации векторов.
Аналогично
$ hat(Vect(y)) = hat(beta)_0 Vect(1) + hat(beta)_1 Vect(x) $

#grid(columns: 2, column-gutter: 0.5em)[
  Так как в задаче поиска МНК оценок мы минимизируем суммарное квадратичное отклонение $y_i - hat(y)_i$, равное квадрату длины вектора $Vect(y) - hat(Vect(y))$, причем вектор $hat(Vect(y))$ лежит в плоскости, определяемой векторами $Vect(1)$ и $Vect(x)$, то можно увидеть, что задача эквивалентна задаче поиска ортогональной проекции $Vect(y)$ на данную плоскость.
  При этом $hat(Vect(y))$ --- проекция $Vect(y)$, $hat(Vect(u)) = Vect(y) - hat(Vect(y))$.
  Нетрудно показать, что $macron(y)$ является длиной проекции векторов $Vect(y)$ и $hat(Vect(y))$ на $Vect(1)$.
][
  #align(center)[#image("fig/OLS.pdf")]
]

Отобразив все приведенные выше факты на графике, мы получим то, что $italic("TSS") = italic("ESS") + italic("RSS")$ является теоремой Пифагора.

= Критерии качества оценивания

Для сравнения качества приближения реальных данных моделью введем два показателя: коэффициент детерминации $R^2$ и стандартную ошибку регрессии $italic("SER")$.

*Коэффициент детерминации* --- доля выборочной дисперсии $Y$, которая объясняется оцененной моделью.
$ R^2 = frac(italic("ESS"), italic("TSS")) = 1 - frac(italic("RSS"), italic("TSS")) $
Несложно показать, что
- $0 lt.eq R^2 lt.eq 1$;
- $R^2 = hat(beta)_1^2 s_X^2 / s_Y^2 = frac(s_(Y X)^2, s_X^2 s_Y^2) = r_(Y X)^2$;
- $R^2 = 0 arrow.l.r.double hat(beta)_1 = 0$.

*Стандартная ошибка регрессии* --- оценка стандартного отклонения случайных ошибок регрессии $u_i$.
$ italic("SER") = s_(hat(u)_i) = sqrt(frac(1, n - 2) sum_(i = 1)^n hat(u)_i^2) $

= Тестирование гипотез

После оценивания модели и проверки того, что получаемые оценки *BLUE*, нам бы хотелось бы протестировать некоторые гипотезы относительно *истинных* значений параметров генеральной совокупности.

Допустим, мы хотим проверить утверждение, что истинное значение $beta_1$ равно некоторому значению $beta_(1,0)$;
то есть протестировать нулевую гипотезу $ H_0 : beta_1 = beta_(1, 0) $ против
- двухсторонней альтернативы $H_a : beta_1 eq.not beta_(1, 0)$ или
- левосторонней альтернативы $H_a : beta_1 < beta_(1, 0)$ или
- правосторонней альтернативы $H_a : beta_1 > beta_(1, 0)$.

Тестирование данной гипотезы проводится при помощи t-теста, о котором вы можете почитать в ликбезе про тестирование гипотез.
Для проведения теста нам нужно получить оценки стандартных ошибок оценок коэффициентов.
Их можно получить при помощи следующих формул:
$
  hat(italic("Var"))(hat(beta)_1 vert X_1\,dots.h) & = underbrace(frac(s_u^2, sum_(i = 1)^n (X_i - macron(X))^2), sigma_u^2 = upright("const")) = underbrace(1 / n times frac(frac(1, n - 2) sum_(i = 1)^n (X_i - macron(X))^2 hat(u)_i^2, [1 / n sum_(i = 1)^n (X_i - macron(X))^2]^2), sigma_u^2 eq.not upright("const"))
$
#proof[$
    hat(beta)_1 - beta_1 & = sum_i hat(a)_i u_i\
    Var [hat(beta)_1 vert X_1, dots.h] & = Ex [(hat(beta)_1 - beta_1)^2 vert X_1, dots.h] = E [sum_i hat(a)_i u_i] =\
    & = sigma_u^2 sum_i hat(a)_i^2 = sigma_u^2 sum frac((X_i - macron(X))^2, sum (X_i - macron(X))^2) = frac(sigma_u^2, sum (X_i - macron(X))^2)
  $ #qedhere]

$
  hat(italic("Var"))(hat(beta)_0 vert X_1, dots.h) & = underbrace(frac((1 / n sum_(i = 1)^n X_i^2) s_u^2, sum_(i = 1)^n (X_i - macron(X))^2), sigma_u^2 = upright("const")) = underbrace(1 / n times frac(frac(1, n - 2) sum_(i = 1)^n hat(H)_i^2 hat(u)_i^2, [1 / n sum_(i = 1)^n hat(H)_i^2]^2), sigma_u^2 eq.not upright("const"))\
  hat(H)_i & = 1 - frac(macron(X), 1 / n sum_(i = 1)^n X_i^2) X_i
$
#proof[$
    Var [hat(beta)_0 vert X_1 \, dots.h] & = Var [macron(Y) - hat(beta)_1 macron(X) vert X_1, dots.h] =\
    & = Var_X macron(Y) + macron(X) Var_X hat(beta)_1 + 2 Cov_X (macron(Y) \, macron(beta)_1) = bigast \
    Var_X macron(Y) & = Var_X [1 / n sum (beta_0 + beta_1 X_i + u_i)] = 1 / n^2 Var_X (sum u_i) = 1 / n sigma_u^2\
    Cov_X (macron(Y) \, macron(beta)_1) & = 1 / n Cov_X [sum Y_i \, sum hat(a)_i Y_i] =\
    & = 1 / n sum sum hat(a)_i underbrace(Cov_X [Y_i, Y_j], =0\, med i eq.not j) = sigma_Y^2 sum hat(a)_i = 0\
    bigast & = sigma_u^2 [1 / n + frac(macron(X), sum (X_i - macron(X))^2)] = sigma_u^2 frac(sum (X_i - macron(X))^2 + n macron(X), n sum (X_i - macron(X))^2) =\
    & = sigma_u^2 frac(sum X_i^2, n sum (X_i - macron(X))^2)
  $ #qedhere]

== Доверительный интервал для $Delta Y$

В соответствии с интерпретацией коэффициента наклона регрессии, при изменении $X$ на $Delta X$ $Y$ изменяется на $beta_1 Delta X$.
Доверительный интервал для изменения $Delta Y$ равен
$ [(hat(beta)_1 - SE hat(beta)_1 times t^"crit") Delta X; (hat(beta)_1 + SE hat(beta)_1 times t^"crit") Delta X] $

= Специальные виды регрессий

Рассмотрим некоторые виды моделей, имеющих особенности в оценивании и интерпретации

== Регрессия без свободного члена

Если в регрессии нет свободного члена, то есть она имеет вид
$ Y_i = beta_1 X_i + u_i $ то:
#grid(columns: (2fr, 1fr))[
  - $hat(beta)_1 eq.not s_(Y X) / s_X^2$;
  - $sum_(i = 1)^n hat(u)_i eq.not 0$;
  - $(macron(Y), macron(X)) in.not hat(Y) = hat(beta)_1 X$;
  - $italic("TSS") eq.not italic("ESS") + italic("RSS")$, то есть $R^2$ не работает!
    В этом случае используем нецентрированный $R^2$ $ R_u^2 = 1 - frac(italic("RSS"), sum_(i = 1)^n Y_i^2) $
][
  #align(center)[#image(height: 30%, "fig/OLS_nocons.pdf")]
]

== Регрессия с бинарной объясняющей переменной

Рассмотрим вопрос интерпретации коэффициента при бинарной переменной, то есть принимающей два значения: 0 и 1.
Для этого предположим, что в нашей выборке размера $n$ переменная $X_i in {0, 1}$, причем число единиц равно $n_1$.

Рассмотрим модель следующего вида
$ Y_i = beta_0 + beta_1 X_i + u_i $

Нетрудно увидеть, что $sum_(i=1)^n X_i = n_1$, и $macron(X) = n_1 / n$.
Также очевидно, что $sum_(i=1)^n X^2_i = n_1$.
Тогда
$
  sum_(i = 1)^n (X_i - macron(X))^2 & = sum_(i = 1)^n X_i^2 - n macron(X) =\
  & = n_1 - n n_1^2 / n^2 = frac(n_1 (n - n_1), n) = frac(n_0 n_1, n) \
  sum_(i = 1)^n (X_i - macron(X))^2 & = sum_(i = 1)^n X_i^2 - n macron(X) =\
  & = n_1 - n n_1^2 / n^2 = frac(n_1 (n - n_1), n) = frac(n_0 n_1, n) \
  sum_(i = 1)^n (X - macron(X))(Y_i - macron(Y)) & = sum_(i = 1)^n X_i (Y_i - macron(Y)) =\
  & = sum_(i = 1)^n X_i Y_i - macron(Y) sum_(i = 1)^n X_i = sum_(i = 1)^n X_i Y_i - n_1 macron(Y) =\
  & = n_1 macron(Y)_1 - n_1 macron(Y)
$
где $macron(Y)_1$ --- выборочное среднее $Y_i$, соответствующих $X_i = 1$.

Тогда имеем
$
  hat(beta)_1 & = frac(n_1 (macron(Y)_1 - macron(Y)), n_0 n_1\/n) = n / n_0 (macron(Y)_1 - macron(Y)) =\
  & = n / n_0 autolr(macron(Y)_1 underbrace(- n_1 / n macron(Y)_1 - n_0 / n macron(Y)_0, macron(Y) = 1 / n (sum_(X_i = 1) Y_i + sum_(X_i = 0) Y_i))) = macron(Y)_1 - macron(Y)_0 \
  hat(beta)_0 & = macron(Y) - hat(beta)_1 macron(X) = macron(Y) - (macron(Y)_1 - macron(Y)_0) n_1 / n =\
  & = n_1 / n macron(Y)_1 + n_0 / n macron(Y)_0 - n_1 / n macron(Y)_1 + n_1 / n macron(Y)_0 = macron(Y)_0
$

То есть в такой модели оценка свободного члена равна выборочному среднему $Y_i$ для $X_i = 0$;
это значение можно интерпретировать как "базовый" случай.
Оценка коэффициента при $X_i$ равна разнице выборочных средних $Y_i$ для "небазового" и "базового" уровней.

== Регрессии с логарифмированными переменными

Часто приходится иметь дело со следующей картиной.

#align(center)[#image(height: 25%, "fig/nonlinear.pdf")]

Понятно, что "нулевая" предпосылка МНК о линейности модели здесь не работает.
Часто спасти ситуацию помогает логарифмирование регрессора и/или регрессанда, приводящее нас в изменной выборке, на которой предположение о линейной связи выполняется.

Вспомним, что для малых $x$ выполняется примерное равенство
$ ln (1 + x) approx x => ln ((x + Delta x)/ x) approx (Delta x) / x $
С учетом этого меняются интерпретации коэффициентов регрессии:

- Линейно-логарифмическая модель
  $ Y_i = beta_0 + beta_1 ln X_i + u_i $
  1%-ное изменение $X$ приводит к изменению $Y$ на $0.01 beta_1$.

- Логарифмически-линейная модель
  $ ln Y_i = beta_0 + beta_1 X_i + u_i $
  Единичное изменение $X$ (то есть $Delta X = 1$) приводит к изменению $Y$ на $100 beta_1%$.

  Важно отметить, что в данном случае нужно быть осторожным с интерпретацией коэффициентов, больших, чем $0.1$.
  Действительно

  #align(center)[#table(
    columns: 2,
    table.hline(),
    table.header[$beta_1$][$exp {beta_1}$],
    table.hline(),
    [0.01], [1.01005],
    [0.02], [1.02020],
    [0.03], [1.03046],
    [0.05], [1.05127],
    [0.09], [1.09417],
    table.hline(),
    [0.1], [1.10520],
    [0.2], [1.22140],
    [0.3], [1.35000],
    [0.5], [1.65000],
    table.hline(),
  )]

- Модель, линейная в логарифмах
  $ ln Y_i = beta_0 + beta_1 ln X_i + u_i $
  1\%-ное изменение $X$ приводит к изменению $Y$ на $beta_1%$.

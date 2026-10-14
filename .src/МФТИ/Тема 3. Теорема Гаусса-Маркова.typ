#import "@preview/theorion:0.6.0": *
#import "@preview/tschich:0.2.0": *
#import "@local/neat-document:0.1.0": *
#import "@local/trinkets:0.1.0": *
#set document(title: [Теорема Гаусса-Маркова])
#show: neat-document.with(margin: tschich-var(210mm, 297mm, 1 / 12))
#show: show-theorion

Ранее мы рассмотрели три предпосылки МНК для парной линейной регрессии.
Напомним:
+ $Ex (u_i mid(bar.v) X_i) = 0$;
+ $(Y_i, X_i) tilde.op "iid"$;
+ $Y_i$ и $X_i$ имеют конечные четвертые моменты --- большие выбросы маловероятны.

При выполнении этих предпосылок МНК оценки коэффициентов регрессии на больших выборках имеют асимптотическое нормальное распределение.
На маленьких выборках для нормального распределения оценок коэффициентов требуется нормальность распределения случайных ошибок.

= Теорема Гаусса-Маркова

Следующая теорема дает нам условия, при выполнении которых полученные оценки коэффициентов парной линейной регрессии будут "хорошими".

== Случай парной регрессии

#theorem[Гаусса-Маркова][
  Рассмотрим уравнение $Y_i = beta_0 + beta_1 X_i + u_i$.
  Если для всех $i, j = 1, dots, n$ выполняются *условия Гаусса-Маркова*:
  + $Ex [u_i mid(bar.v) X_1, dots, X_n] = 0$,
  + $Var [u_i mid(bar.v) X_1, dots, X_n] = sigma_u^2 < oo$,
  + $Cov [u_i, u_j mid(bar.v) X_1, dots, X_n] = 0, quad i eq.not j$,\
  то МНК оценка $hat(beta)_1$ вляется наилучшей (эффективной) линейной несмещенной оценкой (*BLUE*).
]

Перечисленные выше условия Гаусса-Маркова вытекают из предположений МНК:
+ Первое условие Гаусса-Маркова следует из первых двух предположений МНК.
  Действительно, второе условие влечет независимость $u_i$ от $X_1,dots, X_n$, поскольку $ u_i = Y_i - beta_0 - beta_1 X_i $ и все наблюдения независимы.
  Поэтому
  $ Ex [u_i mid(bar.v) X_i] = Ex [u_i mid(bar.v) X_1, dots, X_n] = 0 $
+ Второе условие Гаусса-Маркова следует из второго предположения МНК и _дополнительного_ предположения о *гомоскедастичности* ошибок.

  Случайная ошибка $u_i$ называется *гомоскедастичной*, если $forall i, med Var[u_i mid(bar.v) X_i] = "const"$.
  В противном случае ошибка называется *гетероскедастичной*.
+ Третье условие Гаусса-Маркова следует из второго предположения МНК.

#proof[
  Для доказательства теоремы необходимо определить, какие оценки мы будем называть линейными.

  / Линейная оценка: оценка, имеющая вид $sum_(i=1)^n a_i Y_i$, где $a_i$ --- коэффициенты, *не зависящие* от $Y_i$.

  Вспомним, что
  $ hat(beta)_1 = frac(sum_(i = 1)^n (Y_i - macron(Y))(X_i - macron(X)), sum_(i = 1)^n (X_i - macron(X))^2) = ast.op $
  Раскроем первую скобку в числителе и получим, что
  $
    ast.op = frac(sum_(i = 1)^n (X_i - macron(X))Y_i - macron(Y) sum_(i = 1)^n (X_i - macron(X)), sum_(i = 1)^n (X_i - macron(X))^2) = ast.op
  $
  Так как $sum_(i=1)^n (X_i - macron(X)) = 0$, то
  $
    ast.op = frac(sum_(i = 1)^n (X_i - macron(X))Y_i, sum_(i = 1)^n (X_i - macron(X))^2) = sum_(i = 1)^n frac(X_i - macron(X), sum_(i = 1)^n (X_i - macron(X))^2) Y_i = sum_(i = 1)^n hat(a)_i Y_i
  $
  где $hat(a)_i = frac(X_i - macron(X), sum_(i = 1)^n (X_i - macron(X))^2)$.
  Тем самым мы показали, что $hat(beta)_1$ --- линейная оценка.

  Докажем условную несмещенность $hat(beta)_1$.
  Требуется показать, что $Ex (hat(beta)_1 mid(bar.v) X_1, dots.h, X_n) = beta_1$.
  $
    Ex [hat(beta)_1 mid(bar.v) X_1, dots.h, X_n) & = Ex [beta_1 + frac(sum_(i = 1)^n (X_i - macron(X))u_i, sum_(i = 1)^n (X_i - macron(X))^2) mid(bar.v) X_1, dots.h, X_n] =\
    & = beta_1 + Ex [frac(sum_(i = 1)^n (X_i - macron(X)) u_i, sum_(i = 1)^n (X_i - macron(X))^2) mid(bar.v) X_1, dots.h, X_n] =\
    & = beta_1 + frac(sum_(i = 1)^n (X_i - macron(X)) Ex [u_i mid(bar.v) X_1, dots.h, X_n], sum_(i = 1)^n (X_i - macron(X))^2) =\
    & = beta_1 quad "по 1 условию Г.-М."
  $

  Покажем, что МНК оценка $hat(beta)_1$ является самой эффективной среди всех линейных условно насмещенных оценок.
  Рассмотрим $tilde(beta)_1 = sum_(i=1)^n tilde(a)_i Y_i$ --- *произвольную* линейную условно несмещенную оценку $beta_1$.
  Имеем
  $
    tilde(beta)_1 & = beta_0 sum_(i = 1)^n tilde(a)_i + beta_1 sum_(i = 1)^n tilde(a)_i X_i + sum_(i = 1)^n tilde(a)_i u_i arrow.r.double\
    Ex [tilde(beta)_1 mid(bar.v) X_1, dots.h] & = beta_0 underbrace(sum_(i = 1)^n tilde(a)_i, = 0) + beta_1 underbrace(sum_(i = 1)^n tilde(a)_i X_i, = 1) = beta_1
  $
  В общем случае, для выполнения последнего равенства необходимо выполнение $sum_(i=1)^n tilde(a)_i = 0$ и $sum_(i=1)^n tilde(a)_i X_i = 1$.
  В силу произвольности выбранной оценки данные равенства справедливы для любой линейной условно несмещенной оценки, в том числе и МНК оценки $hat(beta)_1$.

  Найдем дисперсию $tilde(beta)_1$:
  $
    tilde(beta)_1 - beta_1 & = sum_(i = 1)^n tilde(a)_i u_i\
    Var [tilde(beta)_1 mid(bar.v) X_1, dots.h] & = Var [sum_(i = 1)^n tilde(a)_i u_i mid(bar.v) X_1, dots.h] =\
    & = sum_(i = 1)^n sum_(j = 1)^n tilde(a)_i tilde(a)_j Cov [u_i, u_j mid(bar.v) X_1, dots.h] =\
    & = sigma_u^2 sum_(i = 1)^n tilde(a)_i^2 quad "по 2 и 3 условиям Г.-М."
  $

  Предположим, что $tilde(a)_i = hat(a)_i + d_i$, где $hat(a)_i$ --- линейные веса МНК оценки $hat(beta)_1 = sum_(i=1)^n hat(a)_i Y_i$.
  Тогда
  $
    sum_(i = 1)^n tilde(a)_i^2 & = sum_(i = 1)^n hat(a)_i + sum_(i = 1)^n d_i^2 + 2 sum_(i = 1)^n hat(a)_i d_i\
    sum_(i = 1)^n hat(a)_i d_i & = sum_(i = 1)^n frac((X_i - macron(X)) d_i, sum_(i = 1)^n (X_i - macron(X))^2) =\
    & = frac(sum_(i = 1)^n X_i d_i - macron(X) sum_(i = 1)^n d_i, sum_(i = 1)^n (X_i - macron(X))^2) =\
    & = frac(sum_(i = 1)^n tilde(a)_i X_i - sum_(i = 1)^n hat(a)_i X_i - macron(X) (sum_(i = 1)^n tilde(a)_i - sum_(i = 1)^n hat(a)_i), sum_(i = 1)^n (X_i - macron(X))^2) = 0\
    Var [tilde(beta)_1 mid(bar.v) X_1, dots.h] & = sigma_u^2 sum_(i = 1)^n tilde(a)_i^2 = sigma_u^2 sum_(i = 1)^n hat(a)_i^2 + sigma_u^2 sum_(i = 1)^n d_i =\
    & = Var[hat(beta)_1 mid(bar.v) X_1, dots.h] + sigma_u^2 sum_(i = 1)^n d_i > Var[hat(beta)_1 mid(bar.v) X_1, dots.h] quad upright("если") thick tilde(beta)_1 eq.not hat(beta)_1
  $
  Последнее неравенство следует из того, что если $tilde(beta)_1 eq.not hat(beta)_1$ и $tilde(a)_i = hat(a)_i + d_i$, то $exists i : d_i eq.not 0 arrow.l.r.double sum_(i = 1)^n d_i^2 > 0$.
]

== Случай множественной регрессии

#theorem[Гаусса-Маркова][
  Если для модели
  #math.equation(block: true, numbering: "(1)", [$Vect(Y) = Vect(X) Vect(beta) + Vect(u)$]) <eq0>
  выполняются следующие условия
  + $Ex [Vect(u) mid(bar.v) Vect(X)] = 0$;
  + $Sigma_u := Ex[Vect(u)Vect(u)' mid(bar.v) Vect(X)] = sigma^2_u Vect(I)_n$;
  + $rank (Vect(X)'Vect(X)) = k + 1$.
  то МНК оценка $hat(Vect(beta))$ вляется наилучшей (эффективной) линейной несмещенной оценкой (*BLUE*).
]

#proof[
  Линейность $hat(Vect(beta))$ очевидна.
  Действительно, так как
  $ hat(Vect(beta)) = (Vect(X)'Vect(X))^(-1) Vect(X)' Vect(Y) $
  то несложно увидеть, что каждый компонент $hat(Vect(beta))$ представляет собой линейную комбинацию $Y_i$ с весами, определяемыми матрицей $Vect(X)$.
  Также отметим, что
  $
    hat(Vect(beta)) = (Vect(X)'Vect(X))^(-1) Vect(X)' (Vect(X)Vect(beta) + Vect(u)) = Vect(beta) + (Vect(X)'Vect(X))^(-1) Vect(X)' Vect(u)
  $

  Несмещенность также несложно показать:
  $
    Ex[hat(Vect(beta)) mid(bar.v) Vect(X)] & = Ex [(Vect(X)'Vect(X))^(-1) Vect(X)' Vect(Y) mid(bar.v) Vect(X)] = \
    & = Ex [(Vect(X)'Vect(X))^(-1) Vect(X)' (Vect(X)Vect(beta) + Vect(u)) mid(bar.v) Vect(X)] =\
    & = Vect(beta) + (Vect(X)'Vect(X))^(-1) Vect(X)' underbrace(Ex[Vect(u) mid(bar.v) Vect(X)], =0) = Vect(beta) quad "по 1-му условию"
  $

  Перейдем к эффективности#footnote[См. Магнус, Катышев, Пересецкий "Эконометрика. Начальный курс".].
  Пусть $tilde(Vect(beta))$ --- *произвольная* линейная условно насмещенная оценка $Vect(beta)$.
  Ее можно представить в виде
  $ tilde(Vect(beta)) = (Vect(B) + Vect(C))Vect(Y) $
  где $Vect(B) = (Vect(X)'Vect(X))^(-1) Vect(X)'$.

  Из несмещенности $tilde(Vect(beta))$
  $
    Vect(beta) = Ex[tilde(Vect(beta)) mid(bar.v) Vect(X)] = (Vect(B) + Vect(C)) Ex[Vect(Y) mid(bar.v) Vect(X)] = (Vect(B) + Vect(C)) Vect(X) Vect(beta) = (Vect(I)_(k+1) + Vect(C)Vect(X)) Vect(beta)
  $
  следует, что $Vect(C)Vect(X) = Vect(0)$.

  Также отметим, что
  $
    tilde(Vect(beta)) & = (Vect(B) + Vect(C)) Vect(Y) = (Vect(B) + Vect(C)) (Vect(X)Vect(beta) + Vect(u)) =\
    & = underbrace(Vect(B)Vect(X), =Vect(I)_(k+1)) Vect(beta) + underbrace(Vect(C)Vect(X), =Vect(0)) Vect(beta) + (Vect(B) + Vect(C)) Vect(u) =\
    & = Vect(beta) + (Vect(B) + Vect(C)) Vect(u)
  $

  Ковариационная матрица $hat(Vect(beta))$ равна
  $
    Sigma_hat(Vect(beta)) = Ex[(hat(Vect(beta)) - Vect(beta))(hat(Vect(beta)) - Vect(beta))' mid(bar.v) Vect(X)] & = Ex[(Vect(X)'Vect(X))^(-1) Vect(X)'Vect(u)Vect(u)'Vect(X)(Vect(X)'Vect(X))^(-1) mid(bar.v) Vect(X)] =\
    & = (Vect(X)'Vect(X))^(-1) Vect(X)' underbrace(Ex[Vect(u)Vect(u)' mid(bar.v) Vect(X)], Sigma_u = sigma^2_u Vect(I)_n) Vect(X) (Vect(X)'Vect(X))^(-1) =\
    & = sigma^2_u (Vect(X)'Vect(X))^(-1) Vect(X)' Vect(I)_n Vect(X) (Vect(X)'Vect(X))^(-1) =\
    & = sigma^2_u (Vect(X)'Vect(X))^(-1)
  $

  Ковариационная матрица $tilde(Vect(beta))$ равна
  $
    Sigma_tilde(Vect(beta)) = Ex[(tilde(Vect(beta)) - Vect(beta))(tilde(Vect(beta)) - Vect(beta))' mid(bar.v) Vect(X)] & = (Vect(B) - Vect(C)) underbrace(Ex[Vect(u)Vect(u)' mid(bar.v) Vect(X)], Sigma_u = sigma^2_u Vect(I)_n) (Vect(B)' - Vect(C)') =\
    & = sigma^2_u autolr(Vect(B)Vect(B)' + underbrace(Vect(B)Vect(C)' + Vect(C)Vect(B)', =Vect(0) arrow.l Vect(C)Vect(X) = Vect(0)) + Vect(C)Vect(C)') = \
    & = Sigma_hat(Vect(beta)) + sigma^2_u Vect(C)Vect(C)'
  $

  Матрица $Vect(C)Vect(C)'$ --- симметрична и положительно полуопределена (или неотрицательно определена), то есть $forall Vect(x)$
  $ Vect(x)'Vect(C)Vect(C)'Vect(x) >= 0 $
  Это означает, что диагональные элементы матрицы неотрицательны.
  А так как на главных диагоналях ковариационных матриц располагаются дисперсии оценок отдельных коэффициентов, то это означает, что
  $ forall i in {0, 1,dots, k} med Var[tilde(beta)_i] gt.eq Var[hat(beta)_i] $
]

= Нарушение первого условия

Первое условие Гаусса-маркова нарушается, если
$ Ex [u_i mid(bar.v) X_1,dots, X_k] eq.not 0 $
или, эквивалентно в матричной форме,
$ Ex [Vect(u)' Vect(X)] eq.not Vect(0) $
То есть, в этом случае нарушается условие отсутствия корреляции между случайной ошибкой и регрессорами модели.

На лекции про модель множественной линейной регрессии было показано, что в этом случае оценки получаются не только смещенными, но и не состоятельными.
Одним из источников нарушения первого УГМ является пропуск значимой объясняющей переменной, коррелированной с включенным регрессором.
На самом деле, это не единственный источник, но подробнее о них мы поговорим на последнем занятии.

= Нарушение второго условия

Напомним, что второе условие Гаусса-Маркова требует
$ Var [u_i mid(bar.v) X_1, dots.h, X_k] = sigma_u^2 < oo $
или в матричной форме (если объединить с третьим условием)
$ Sigma_u = sigma_u^2 Vect(I)_(n times n) $

При невыполнении этого условия дисперсии случайных ошибок для различных наблюдений будут неодинаковыми.
Это называется *гетероскедастичностью*.
#image("fig/fig_6_1.png")

Это приведет к тому, что оценки дисперсий оценок будут смещены.
При этом сами оценки останутся несмещенными и состоятельными.
То есть
- оценки станут неэффективными,
- статистические тесты будут смещенными.

== Взвешенный МНК

"Наилучший" с точки зрения эффективности оценок способ --- взвешенный метод наименьших квадратов (ВМНК).

Рассмотрим парную модель
#math.equation(block: true, numbering: "(1)", [$Y_i = beta_0 + beta_1 X_i + u_i$]) <eq1>
Пусть
$ Var[u_i mid(bar.v) X_i] = lambda h(X_i) $
где $lambda$ --- некоторая константа, $h(X_i)$ --- _известная_ функция.
Рассмотрим преобразованные переменные
$
  tilde(Y)_i = Y_i / sqrt(h(X_i)), quad tilde(bb(1))_i = 1 / sqrt(h(X_i)), quad tilde(X)_i = X_i / sqrt(h(X_i)), quad tilde(u)_i = u_i / sqrt(h(X_i))
$
и подставим их в уравнение @eq1 вместо исходных.
Тогда модель примет вид
#math.equation(
  block: true,
  numbering: "(1)",
  [$tilde(Y)_i = beta_0 tilde(bb(1))_i + beta_1 tilde(X)_i + tilde(u)_i$],
) <eq2>
$hat(beta)_1$ --- *ВМНК* оценка коэффициента наклона в уравнении @eq1, и она является BLUE.

Иначе ВМНК можно выразить в терминах матриц.
Так как $Sigma_u = op("diag"){sigma^2_1, dots, sigma^2_N}$, то обратная к ней матрица будет иметь вид
$ Sigma^(-1)_u = op("diag"){1/sigma^2_1, dots, 1/sigma^2_N} $
Данную матрицу можно представить в виде $Sigma^(-1)_u = Vect(L)'Vect(L)$, где $Vect(L) = op("diag"){1/sigma_1, dots, 1/sigma_N}$.

Если домножить уравнение @eq0 на матрицу $Vect(L)$ слева, то получим оценку
$ hat(Vect(beta))^"WLS" = (Vect(X)' Sigma^(-1)_u Vect(X))^(-1) Vect(X)' Sigma^(-1)_u Vect(y) $
которая будет *BLUE*.

Главная проблема состоит в том, что точный вид функции $h(X_i)$ нам неизвестен; по этой причине полученная выше оценка называется _недоступной_.
Для получения _доступной_ оценки нужно
+ Оцените регрессию @eq1 при помощи МНК и вычислите остатки $hat(u)_i$.
+ Оцените _предполагаемую_ модель функции условной дисперсии $Var[u_i mid(bar.v) X_i]$: оцените регрессию $hat(u)_i^2$ на набор переменных, от которых она зависит.
+ Используйте оцененную модель для расчета предсказанных значений функции условной дисперсии $hat(Var)[u_i mid(bar.v) X_i]$.
+ Умножьте зависимую переменную и регрессоры (в том числе и свободный член) на величину, обратную квадратному корню из оцененной функции условной дисперсии $hat(Var)[u_i mid(bar.v) X_i]$.
+ Оцените коэффициенты взвешенной регрессии @eq2 при помощи МНК; полученные оценки являются ВМНК–оценками @eq1.

Главная проблема ВМНК: как выбрать спецификацию $Var[u_i mid(bar.v) X_i]$?
Дело в том, что неверный выбор спецификации сохраняет и, иногда, даже усубляет неэффективность оценок.

== Устойчивые к гетероскедастичности ошибки Уайта

Устойчивые к гетероскедастичности стандартные ошибки Эйкера-Хьюбера-Уайта (состоятельные при выполнении трех предположений МНК):
$
  hat(sigma)_(hat(beta)_1)^2 & = 1 / n frac(frac(1, n - 2) sum_(i = 1)^n (X_i - macron(X))^2 hat(u)_i^2, [1 / n sum_(i = 1)^n (X_i - macron(X))^2]^2)\
  hat(sigma)_(hat(beta)_0)^2 & = 1 / n frac(frac(1, n - 2) sum_(i = 1)^n hat(H)_i^2 hat(u)_i^2, [1 / n sum_(i = 1)^n hat(H)_i^2]^2)
$
где $hat(H)_i = 1 - frac(macron(X), 1 / n sum_(i = 1)^n X_i^2) X_i$.

Для множественной регрессии формула имеет вид
$
  hat(Sigma)_hat(Vect(beta)) = (Vect(X)'Vect(X))^(-1)(Vect(X)'op("diag"){hat(u)^2_1, dots, hat(u)^2_n} Vect(X))(Vect(X)'Vect(X))^(-1)
$

Данные оценки _все еще_ неэффективны, но в меньшей степени, чем при их неиспользовании.
На практике этот метод спользуется наиболее широко, зачастую "по умолчанию".

== Тестирование гетероскедастичности

=== Тест Голдфельдта-Квандта

Применяется в предположении, что дисперсия зависит от одного из регрессоров.

Для проведения теста надо отсортировать выборку по предположительному возрастанию дисперсии ошибок.
После чего надо отбросить $r$ наблюдений из середины выборки и оценить исходную модель на каждой из оставшихся подвыборок отдельно, получив оценки $italic("RSS")$.
После чего надо протестировать нулевую гипотезу
$ H_0: italic("RSS")_1 = italic("RSS")_2 $
против альтернативы
$ H_a: italic("RSS")_1 < italic("RSS")_2 $
Для этого воспользуемся статистикой
$ F = frac(italic("RSS")_2, italic("RSS")_1) tilde.op F_(1 - alpha) (frac(n - r, 2) - 2 \; frac(n - r, 2) - 2) $

=== Тесты Бройша-Пагана и Уайта

Сущность тестов состоит в следующем.
Сначала оценивается базовая модель
$ Vect(Y) = Vect(X) Vect(beta) + Vect(u) $
из которой находятся остатки $hat(u)_i$.

Затем для теста Бройша-Пагана оценивается вспомогательная регрессия
$ hat(Vect(u))^2 = gamma_0 + sum_i gamma_i Vect(z)_i + Vect(nu) $
где $Vect(Z)$ --- матрица значений переменных, которые, как мы считаем, влияют на дисперсию остатков. Для теста Уайта эта регрессия имеет вид
$ hat(Vect(u))^2 = gamma_0 + sum_i gamma_i Vect(z)_i + sum_i sum_(j>=i) gamma_(i j) Vect(z)_i Vect(z)_j + Vect(nu) $

После необходимо протестировать гипотезу об одновременной незначимости коэффициентов $Vect(gamma)$ при помощи статистик
$
  n R^2 & tilde.op chi^2(p) \
      F & tilde.op F(n - 1, n - k - 1)
$
где $R^2$ --- коэффициент детерминации во вспомогательной регрессии, $p$ --- число регрессоров во вспомогательной регрессии.

= Нарушение третьего условия

Нарушение третьего условия Гаусса-Маркова означает, что
$ exists i, j: Ex[u_i u_j mid(bar.v) X_1, dots, X_k] eq.not 0 $
или, в матричной форме, наличие в ковариационной матрице ошибок ненулевых элементов вне главной диагонали.
Данное явление называется *автокорреляцией*.
Автокорреляция не приводит к смещению или несостоятельности оценок,
но они перестают быть эффективными.

#box(stroke: 2pt + black, inset: 1em)[
  Важное замечание: об автокорреляции имеет смысл говорить только тогда, когда данные имеют *четкий порядок*!
  Для межобъектных данных говорить об автокорреляции нет никакого смысла!
]

== Процедура Кохрейна-Оркатта

Пусть
#math.equation(block: true, numbering: "(1)")[$
      Y_t & = beta_0 bb(1)_t + beta_1 X_t + u_t \
      u_t & = rho u_(t - 1) + nu_t \
  bb(1)_t & = 1 \
        t & = 1\,dots.h\,T \
     nu_t & tilde.op "WN"
$] <eq4>

Рассмотрим
$ rho Y_(t - 1) = beta_0 bb(1)_t rho + beta_1 rho X_(t - 1) + rho u_(t - 1) $

Вычтем выражение из первого выражения второе
#math.equation(
  block: true,
  numbering: "(1)",
)[$underbrace(Y_t - rho Y_(t - 1), tilde(Y_i)_t) = beta_0 underbrace((bb(1)_t - rho bb(1)_(t - 1)), tilde(bb(1))_t) + beta_1 underbrace((X_t - rho X_(t - 1)), tilde(X)_t) + underbrace((u_t - rho u_(t - 1)), nu_t)$] <eq5>

Так как на шаге (3) теряется первое наблюдение, то теоретически возможна ситуация, когда в результате применения процедуры мы получим _менее эффективные_ оценки, чем до нее.
Поэтому применяют поправку Прайса-Уинстена:
$
  tilde(Y)_1 & = sqrt(1 - hat(rho)^2) Y_1 \
     bb(1)_1 & = sqrt(1 - hat(rho)^2) \
  tilde(X)_1 & = sqrt(1 - hat(rho)^2) X_1
$

Поскольку значение $rho$ нам неизвестно, то применяется следующая процедура:
+ Оцениваем модель @eq4 и вычисляем остатки $hat(u)_t$.
+ Оцениваем модель $hat(u)_t = rho hat(u)_t + nu_t$ и получаем $hat(rho)$.
+ Вычисляем $tilde(Y)_t = Y_t - hat(rho) Y_(t-1)$, $tilde(X)_t = X_t - hat(rho) X_(t-1)$, $bb(1)_t = 1 - hat(rho)$.
+ Оцениваем модель @eq5 для $tilde(Y)$, $tilde(X)_t$ и $bb(1)_t$.
+ Возвращаемся к пункту (2) до тех пор, пока аболютная величина разницы между оценками не достигнет заданной точности.

== Обобщенный МНК

При наличии нарушения условия отсутствия автокорреляции ковариационная матрица ошибок перестает быть диагональной. Аналогично подходу ВМНК обратную к ней матрицу можно представиь в виде
$
  Vect(Sigma)_u = Ex[Vect(u)Vect(u)'] = underbrace(frac(sigma_epsilon^2, 1 - rho^2), sigma_u^2) mat(delim: "(", 1, rho, rho^2, dots.h.c, rho^(T - 1); rho, 1, rho, dots.h.c, rho^(T - 2); rho^2, rho, 1, dots.h.c, rho^(T - 3); dots.v, dots.v, dots.v, dots.down, dots.v; rho^(T - 1), rho^(T - 2), rho^(T - 3), dots.h.c, 1) = sigma_u^2 Vect(R)
$
$
  Vect(R)^(-1) = frac(1, 1 - rho^2) mat(delim: "(", 1, - rho, 0, dots.h.c, 0, 0; - rho, 1 + rho^2, - rho, dots.h.c, 0, 0; 0, - rho, 1 + rho^2, dots.h.c, 0, 0; dots.v, dots.v, dots.v, dots.down, dots.v, dots.v; 0, 0, 0, dots.h.c, 1 + rho^2, - rho; 0, 0, 0, dots.h.c, - rho, 1)
$

Матрица $Vect(R)$ единственным образом раскладывается на произведение $Vect(L)'Vect(L)$, где $Vect(L)$ --- нижнетреугольная матрица
$
  Vect(L) = 1 / sqrt(1 - rho^2) mat(delim: "(", sqrt(1 - rho^2), 0, 0, dots.h.c, 0, 0; - rho, 1, 0, dots.h.c, 0, 0; 0, - rho, 1, dots.h.c, 0, 0; dots.v, dots.v, dots.v, dots.down, dots.v, dots.v; 0, 0, 0, dots.h.c, - rho, 1)
$
Из полученной матрицы видно, что поправка Прайса-Уинстена возникает автоматически после того, как уравнение @eq4 домножается слева на $Vect(L)$.

Оценка модифицированного уравнения называется *обобщенной МНК* оценкой
$
  hat(Vect(beta))_"GLS" & = (Vect(X)' Vect(L)' Vect(L) Vect(X))^(- 1) Vect(X)' Vect(L)' Vect(L)Vect(y)=\
  & = (Vect(X)' Vect(R)^(- 1) Vect(X))^(- 1) Vect(X)' Vect(R)^(- 1) Vect(y)=\
  & = (Vect(X)' Vect(Sigma)_(Vect(u))^(- 1) Vect(X))^(- 1) Vect(X)' Vect(Sigma)_(Vect(u))^(- 1) Vect(y)
$
Обратите внимание, что данная оценка аналогина ВМНК оценке; то есть ВМНК оценка является частным случаем ОМНК оценки.

== Тестирование автокорреляции

=== Тест Дарбина-Уотсона

- Применяется только для тестирования автокорреляции первого порядка
  $ u_t = rho u_(t-1) + nu_t $
- В регрессии должен присутствовать свободный член;
- В регрессии не должно быть лаговых значений регрессанда;
- Не должно быть пропусков в переменных.

+ Оцениваем регрессию при помощи МНК и вычисляем остатки $hat(u)_i,med i = 1, dots.h, n$.
+ Вычисляем статистику $ "DW" = frac(sum_(i=2)^T\(hat(u)_t - hat(u)_(t - 1)\)^2, sum_(i = 1)^T hat(u)_t^2) $
+ Отвергаем $H_0: rho = 0$ в пользу $H_a: rho > 0$ на уровне значимости $alpha$ если $"DW" < d_(L alpha)$; не отвергаем, если $"DW" > d_(U alpha)$. Отрезок $d_(L alpha) < "DW" < d_(U alpha)$ представляет собой область неопределенности.

  Для тестирования с использованием альтернативы $H_a: rho < 0$ критические значения надо отразить относительно $2$:
  $
    d^*_(L alpha) & = 4 - d_(L alpha) \
    d^*_(U alpha) & = 4 - d_(U alpha)
  $

=== Тест Бройша-Годфри

Сначала оценивается базовая модель
$ Vect(Y) = Vect(X) Vect(beta) + Vect(u) $
из которой находятся остатки $hat(u)_i$.

Затем оценивается вспомогательная регрессия
$ hat(u)_t = Vect(X)' Vect(beta) + sum_(i = 1)^k gamma_i hat(u)_(t - i) + nu_t $

После необходимо протестировать гипотезу об одновременной незначимости коэффициентов $Vect(gamma)$ при помощи статистик
$
  n R^2 & tilde.op chi^2(k) \
      F & tilde.op F(n - 1, n - k)
$

#import "@preview/theorion:0.6.0": *
#import "@preview/tschich:0.2.0": *
#import "@local/neat-document:0.1.0": *
#import "@local/trinkets:0.1.0": *
#set document(title: [Автокорреляционная функция])
#show: neat-document.with(margin: tschich-var(210mm, 297mm, 1 / 12))
#show: show-theorion

= Aвтокорреляция

*Автокорреляционной функцией* слабостационарного случайного процесса $y_t$ называется функция $rho (s) = rho_s = display(frac(Cov (y_t,y_(t - s)), sqrt(Var y_t Var y_(t - s)))) = display(frac(gamma_s, gamma_0))$.

*Выборочной автокорреляционной функцией* называется функция
$
  r_s = frac(sum^T_(t = s + 1)(y_t - macron(y)) (y_(t - s) - macron(y)), sum^T_(t = 1)(y_t - macron(y))^2)
$
- $r_s$ в случае независимых наблюдений является смещенной оценкой.
- С ростом $s$ объем используемой выборки, а следовательно и число степеней свободы, уменьшается.
  Поэтому рекомендуется ограничить порядок $s$ и не превышать $T/4$.
- Матрица выборочных автокорреляций положительно определена.
- $r_s limits(->)_(T -> infinity) rho_s$.

Автокорреляция, являясь важной характеристикой случайного процесса, обладает определенным недостатком.
Давайте рассмотрим процесс AR(2)
$ y_t = a_1 y_(t - 1) + a_2 y_(t - 2) + epsilon_t $

Если мы найдем значение $rho_2$, то оно, строго говоря, будет показывать нечто среднее между "прямой" корреляцией $y_t$ и $y_(t - 2)$ и "опосредованной" корреляцией, возникающей из-за того, что $y_t$ зависит от величины $y_(t - 1)$, которая, в свою очередь, зависит от $y_(t - 2)$.
#align(center)[#image("fig/PACF.pdf", width: 60%)]

Для того, чтобы обойти этот эффект, используют *частные автокорреляции*.
*Частная автокорреляция* порядка $s$ определяется как значение коэффициента $phi_(s s)$ в следующей модели AR(s):
#eqnum[$ y_t = phi_(s 0) + phi_(s 1)y_(t - 1) + phi_(s 2)y_(t - 2) + dots.h.c+ phi_(s s)y_(t - s) + nu_t $] <eq-6.1>

Обратите внимание, что при поиске значения $phi.alt_(s s)$ мы полностью игнорируем истинный порядок авторегрессии процесса; значение $s$ может принимать произвольные интересующие на значения.

Альтернативным способом расчета частных автокорреляций является использование системы уравнений Юла-Уокера и метода Крамера.
Вернемся к уравнению @eq-6.1.
Построим систему уравнений
#eqnum[$
  rho_j = phi_(s 1)rho_(j - 1) + phi_(s 2)rho_(j - 2) + dots.h.c+ phi_(s s)rho_(j - s),quad j= overline(#[$1, s$])
$] <eq-6.2>

Для того, чтобы найти $phi_(s s)$, воспользуемся правилом Крамера.
Напомним, что значение $phi_(s s)$ равно дроби, в знаменателе которой находится определитель матрицы коэффициентов системы
$
  mat(
    delim: "|", rho_0, rho_1, rho_2, dots.h.c, rho_(s - 2), rho_(s - 1);
    rho_1, rho_0, rho_1, dots.h.c, rho_(s - 3), rho_(s - 2);
    dots.v, dots.v, dots.v, dots.down, dots.v, dots.v;
    rho_(s - 2), rho_(s - 3), rho_(s - 4), dots.h.c, rho_0, rho_(1);
    rho_(s - 1), rho_(s - 2), rho_(s - 3), dots.h.c, rho_1, rho_(0)
  ),
$
а в числителе — определитель матрицы коэффициентов, в которой последний столбец заменен на столбец свободных членов ("правых частей")
$
  mat(
    delim: "|", rho_0, rho_1, rho_2, dots.h.c, rho_(s - 2), rho_(1);
    rho_1, rho_0, rho_1, dots.h.c, rho_(s - 3), rho_(2);
    dots.v, dots.v, dots.v, dots.down, dots.v, dots.v;
    rho_(s - 2), rho_(s - 3), rho_(s - 4), dots.h.c, rho_0, rho_(s - 1);
    rho_(s - 1), rho_(s - 2), rho_(s - 3), dots.h.c, rho_1, rho_(s)
  )
$

*Выборочная частная автокорреляционная функция* $hat(phi)_(s s)$ определяется как коэффициент $phi_(s s)$ в выборочной регрессии
$ (y_t - macron(y)) = phi_(s 1)(y_(t - 1) - macron(y)) + dots.h.c+ phi_(s s)(y_(t - s) - macron(y)) + epsilon_t $

Можно показать, что $hat(phi)_(s s)$ — состоятельная оценка $phi_(s s)$.
Для поиска значения $hat(phi)_(s s)$ можно также воспользоваться правилом Крамера, заменив истинные автокорреляции выборочными.

== AR(1)

Для удобства перейдем к центрированным значениям процесса.
Нетрудно показать, что формула для центрированного авторегрессионного процесса $y_t$ отличается от формулы для исходного процесса отсутствием
свободного члена.

Действительно, пусть мы имеем процесс AR(p)
$ y_t = a_0 + a_1 y_(t - 1) + dots.h.c+ a_p y_(t - p) + epsilon_t $
Нам уже известно, что математическое ожидание этого процесса равно
$ Ex y_t = mu = frac(a_0, 1 - sum^p_(i = 1)a_i) $
Вычтем $mu$ из обеих частей процесса:
$
  macron(y)_t = y_t - mu & = a_0 - frac(a_0, 1 - sum^p_(i = 1)a_i) + a_1 y_(t - 1) + dots.h.c+ a_p y_(t - p) + epsilon_t = \
  & = frac(a_0 - a_0 a_1 - a_0 a_2 - dots.h.c- a_0 a_p - a_0, 1 - sum^p_(i = 1)a_i) + a_1 y_(t - 1) + dots.h.c+ a_p y_(t - p) + epsilon_t = \
  & = frac(- a_0 a_1 - a_0 a_2 - dots.h.c- a_0 a_p, 1 - sum^p_(i = 1)a_i) + a_1 y_(t - 1) + dots.h.c+ a_p y_(t - p) + epsilon_t = bigast
$

Первое слагаемое можно рассматривать как сумму из $p$ отдельных дробей.
Выделим первую из них.
$
  bigast& = frac(- a_0 a_1, 1 - sum^p_(i = 1)a_i) + frac(- a_0 a_2 - dots.h.c- a_0 a_p, 1 - sum^p_(i = 1)a_i) + a_1 y_(t - 1) + dots.h.c+ a_p y_(t - p) + epsilon_t = \
  & = frac(- a_0 a_2 - dots.h.c- a_0 a_p, 1 - sum^p_(i = 1)a_i) + a_1 (y_(t - 1) - frac(a_0, 1 - sum^p_(i = 1)a_i)) + dots.h.c+ a_p y_(t - p) + epsilon_t = \
  & = frac(- a_0 a_2 - dots.h.c- a_0 a_p, 1 - sum^p_(i = 1)a_i) + a_1 underbrace((y_(t - 1) - mu), macron(y)_(t - 1)) + dots.h.c+ a_p y_(t - p) + epsilon_t = bigast
$

Проделав аналогичные действия с остальными слагаемыми, получим
$ bigast= a_1 macron(y)_(t - 1) + dots.h.c+ a_p macron(y)_(t - p) + epsilon_t $

Очевидно, что значения всех автокорреляций для процессов $y_t$ и $macron(y)_t$ совпадают.
Поэтому, для упрощения и для краткости будем считать, что мы уже провели центрирование процесса и отбросили свободный член.

+ Чтобы найти первую автокорреляцию $rho_1$ найдем $gamma_1$, умножив обе части процесса на $y_(t - 1)$ и взяв математическое ожидание.
  Так как мы условились считать, что мы уже провели центрирование, то математическое ожидание процесса равно нулю, и полученное математическое ожидание в очности равно автоковариации.
  $
    gamma_1 = Ex (y_t y_(t - 1)) & = Ex (a_1 y^2_(t - 1) + epsilon_t y_(t - 1)) = \
    & = a_1 underbrace(Ex y^2_(t - 1), = gamma_0) + = underbrace(Ex (epsilon_t y_(t - 1)), 0) = \
    & = a_1 gamma_0 \
    rho_1 = frac(gamma_1, gamma_0) & = a_(1)
  $

  Последнее слагаемое равно нулю, так как мы предполагаем, что прошлое значение процесса $y_(t - 1)$ не испытывают влияния будущей инновации $epsilon_t$.

+ Автокорреляция $rho_2$ ищется аналогичным способом:
  $
    gamma_2 = Ex (y_t y_(t - 2)) & = Ex (a_1 y_(t - 1)y_(t - 2) + epsilon_t y_(t - 2)) = \
    & = a_1 underbrace(Ex (y_(t - 1)y_(t - 2)), = gamma_1) + underbrace(Ex (epsilon_t y_(t - 2)), = 0) = \
    & = a_1 gamma_1 \
    rho_2 = frac(gamma_2, gamma_0) & = a_1 rho_1 = a^2_(1)
  $

  Нетрудно показать, что для процесса AR(1) выполняется
  $ rho_s = a_1 rho_(s - 1) = a^s_1 $

Теперь найдем его частные автокорреляции.

+ Из уравнения
  $ y_t = phi_(01) + phi_(11)y_(t - 1) + nu_t $
  сразу следует, что#footnote[Это справедливо для *любого* слабостационарного процесса.]
  $ phi.alt_(11) = rho_1 $
+ Из уравнения
  $ y_t = phi_(02) + phi_(12)y_(t - 1) + phi_(22)y_(t - 2) + nu_t $
  следует, что
  $
    phi_(22) = frac(
      mat(
        delim: "|", rho_0, rho_(1);
        rho_1, rho_(2)
      ), mat(
        delim: "|", rho_0, rho_(1);
        rho_1, rho_(0)
      )
    ) = frac(mat(delim: "|", 1, a_(1); a_1, a^2_1), mat(delim: "|", 1, a_(1); a_1, 1)) = 0
  $

  Очевидно, что последний столбец определителя в числителе является первым столбцом, умноженным на $a_1$.
  Другими словами, он равен нулю.

  Нетрудно показать, что для $phi_(s s) = 0, forall s > 1$.

== AR(2)

Рассмотрим слабостационарный процесс $y_t = a_0 + a_1 y_(t - 1) + a_2 y_(t - 2) + epsilon_t$.
Как и ранее, упростим себе задачу, центрировав процесс, то есть отбросив свободный член.

+ Найдем первую автокорреляцию:
  $
    gamma_1 = Ex (y_t y_(t - 1)) & = Ex [ a_1 y^2_(t - 1) + a_2 y_(t - 2)y_(t - 1) + epsilon_t y_(t - 1) ] = \
    & = a_1 underbrace(Ex y^2_(t - 1), = gamma_0) + a_2 underbrace(Ex y_(t - 2)y_(t - 1), = gamma_1) + underbrace(Ex epsilon_t y_(t - 1), = 0) = \
    & = a_1 gamma_0 + a_2 gamma_1 \
    rho_1 = frac(gamma_1, gamma_0) & = a_1 + a_2 rho_1 => rho_1 = frac(a_1, 1 - a_2)
  $

+ Найдем вторую автокорреляцию:
  $
    gamma_2 = Ex (y_t y_(t - 2)) & = Ex [ a_1 y_(t - 1)y_(t - 2) + a_2 y^2_(t - 2) + epsilon_t y_(t - 2) ] = \
    & = a_1 underbrace(Ex (y_(t - 1)y_(t - 2)), = gamma_1) + a_2 underbrace(Ex y^2_(t - 2), = gamma_0) + underbrace(Ex (epsilon_t y_(t - 2)), = 0) = \
    & = a_1 gamma_1 + a_2 gamma_0 \
    rho_2 = frac(gamma_2, gamma_0) & = a_1 rho_1 + a_(2)
  $

+ Найдем третью автокорреляцию:
  $
    gamma_3 = Ex (y_t y_(t - 3)) & = Ex [ a_1 y_(t - 1)y_(t - 3) + a_2 y_(t - 2)y_(t - 3) + epsilon_t y_(t - 3) ] = \
    & = a_1 underbrace(Ex (y_(t - 1)y_(t - 3)), = gamma_2) + a_2 underbrace(Ex (y_(t - 2)y_(t - 3)), = gamma_1) + underbrace(Ex (epsilon_t y_(t - 3)), = 0) = \
    & = a_1 gamma_2 + a_2 gamma_1 \
    rho_3 = frac(gamma_3, gamma_0) & = a_1 rho_2 + a_2 rho_(1)
  $

Нетрудно показать, что для процесса AR(2) выполняется
$ rho_s = a_1 rho_(s - 1) + a_2 rho_(s - 2) $

Теперь найдем его частные автокорреляции.
+ Ранее мы уже отмечали, что для *любого* слабостационарного процесса
  $ phi.alt_(11) = rho_1 $

+ Из уравнения
  $ y_t = phi_(02) + phi_(12)y_(t - 1) + phi_(22)y_(t - 2) + nu_t $
  следует, что
  $
    phi_(22) & = frac(
      mat(
        delim: "|", rho_0, rho_(1);
        rho_1, rho_(2)
      ), mat(
        delim: "|", rho_0, rho_(1);
        rho_1, rho_(0)
      )
    ) = frac(rho_2 - rho^2_1, 1 - rho^2_1) = frac(frac(a^2_1 - a^2_2 + a_2, 1 - a_2) - frac(a^2_1, (1 - a_2)^2), 1 - frac(a^2_1, (1 - a_2)^2)) = \
    & = frac(a^3_2 - 2 a^2_2 + a_2 - a^2_1 a_2, (1 - a_2)^2 - a^2_1) = frac(a_(2)[ (1 - 2 a_2 + a^2_2) - a^2_1 ], (1 - a_2)^2 - a^2_1) = a_(2)
  $

+ $
    phi_(33) = frac(
      mat(
        delim: "|", rho_0, rho_1, rho_(1);
        rho_1, rho_0, rho_(2);
        rho_2, rho_1, rho_(3)
      ), mat(
        delim: "|", rho_0, rho_1, rho_(2);
        rho_1, rho_0, rho_(1);
        rho_2, rho_1, rho_(0)
      )
    ) = 0
  $
  Очевидно, что последний столбец определителя в числителе является линейной комбинацией первых двух столбцов с коэффициентами $a_1,a_2$.
  Другими словами, он равен нулю.

Нетрудно показать, что для $phi_(s s) = 0, forall s> 2$.

== MA(1)

Рассмотрим слабостационарный процесс $y_t = a_0 + epsilon_t + beta_1 epsilon_(t - 1)$.
Очевидно, что $mu = a_0$, то есть центрирование также дает нам процесс, чья формула совпадает с формулой исходного процесса, но без свободного члена.

Как и ранее, упростим себе задачу, центрировав процесс, то есть отбросив свободный член.

+ Найдем первую автокорреляцию:
  $
    gamma_0 = Ex y^2_t & = Ex (epsilon_t + beta_1 epsilon_(t - 1))^2 = \
    & = underbrace(Ex epsilon^2_t, = sigma^2) + 2 beta_1 underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta^2_1 underbrace(Ex epsilon^2_(t - 1), = sigma^2) = \
    & = (1 + beta^2_1)sigma^2 \
    gamma_1 = (y_t y_(t - 1)) & = Ex (epsilon_t + beta_1 epsilon_(t - 1)) (epsilon_(t - 1) + beta_1 epsilon_(t - 2)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 1)), = sigma^2) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) = \
    & = beta_1 sigma^2 \
    rho_1 = frac(gamma_1, gamma_0) & = frac(beta_1, 1 + beta^2_1)
  $

  Обратим внимание, что для процесса MA(1) первая автокорреляция не превышает по абсолютному значению $0.5$!

+ Найдем вторую автокорреляцию:
  $
    gamma_2 = Ex (y_t y_(t - 2)) & = Ex (epsilon_t + beta_1 epsilon_(t - 1)) (epsilon_(t - 2) + beta_1 epsilon_(t - 3)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) = \
    & = 0 \
    rho_2 = frac(gamma_2, gamma_0) & = 0
  $

  Все последующие автокорреляции также равны 0.

Теперь найдем его частные автокорреляции.
+ Ранее мы уже отмечали, что для *любого* слабостационарного процесса
  $ phi.alt_(11) = rho_1 $

+ $
    phi_(22) = frac(
      mat(
        delim: "|", rho_0, rho_(1);
        rho_1, 0
      ), mat(
        delim: "|", rho_0, rho_(1);
        rho_1, rho_(0)
      )
    ) = frac(- rho^2_1, 1 - rho^2_1)
  $

+ $
    phi_(33) = frac(
      mat(
        delim: "|", rho_0, rho_1, rho_(1);
        rho_1, rho_0, 0 ;
        0, rho_1, 0
      ), mat(
        delim: "|", rho_0, rho_1, 0 ;
        rho_1, rho_0, rho_(1);
        0, rho_1, rho_(0)
      )
    ) = frac(rho^3_1, 1 - 2 rho^2_1)
  $

+ $
    phi_(44) = frac(
      mat(
        delim: "|", rho_0, rho_1, 0, rho_(1);
        rho_1, rho_0, rho_1, 0 ;
        0, rho_1, rho_0, 0 ;
        0, 0, rho_1, 0
      ), mat(
        delim: "|", rho_0, rho_1, 0, 0 ;
        rho_1, rho_0, rho_1, 0 ;
        0, rho_1, rho_0, rho_(1);
        0, 0, rho_1, rho_(0)
      )
    ) = frac(- rho^4_1, 1 - rho^2_1 + rho^4_1)
  $

Нетрудно заметить, что в знаменателе *всегда* стоит положительное значение, поскольку это определитель корреляционной матрицы, которая положительно определена.
Поэтому знак частной автокорреляции будет зависеть от знака $rho_1$, который зависит от знака $beta_1$:
- Если $beta_(1)> 0$, то
  $
    phi_(s s) = cases(
      > 0\, & s #text[нечетно],
      < 0\, & s #text[четно]
    )
  $
- Если $beta_(1)< 0$, то $phi_(s s)< 0, forall s$.

По абсолютному значению $phi_(s s)$ будет убывать в обоих случаях.

== MA(2)

Рассмотрим слабостационарный процесс $y_t = a_0 + epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2)$.
Как и ранее, упростим себе задачу, центрировав процесс, то есть отбросив свободный член.

+ Найдем первую автокорреляцию:
  $
    gamma_0 = Ex y^2_t & = Ex (epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2))^2 = \
    & = underbrace(Ex epsilon^2_t, sigma^2) + beta^2_1 underbrace(Ex epsilon^2_(t - 1), sigma^2) + beta^2_2 underbrace(Ex epsilon^2_(t - 2), sigma^2) + \
    & wide + 2 beta_1 underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + 2 beta_2 underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + 2 beta_1 beta_2 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) = \
    & = (1 + beta^2_1 + beta^2_2) sigma^(2)
  $
  $
    gamma_1 = Ex (y_t y_(t - 1)) & = Ex (epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2)) (epsilon_(t - 1) + beta_1 epsilon_(t - 2) + beta_2 epsilon_(t - 3)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_2 underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + \
    & wide + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 1)), = sigma^2) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) + beta_1 beta_2 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) + \
    & wide + beta_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 1)), = 0) + beta_2 beta_1 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 2)), = sigma^2) + beta^2_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 3)), = 0) = \
    & = (beta_1 + beta_1 beta_2) sigma^2 \
    rho_1 = frac(gamma_1, gamma_0) & = frac(beta_1 + beta_1 beta_2, 1 + beta^2_1 + beta^2_2)
  $

+ Найдем вторую автокорреляцию:
  $
    gamma_2 = Ex (y_t y_(t - 2)) & = Ex (epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2)) (epsilon_(t - 2) + beta_1 epsilon_(t - 3) + beta_2 epsilon_(t - 4)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + beta_2 underbrace(Ex (epsilon_t epsilon_(t - 4)), = 0) + \
    & wide + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) + beta_1 beta_2 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 4)), = 0) + \
    & wide + beta_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 2)), = sigma^2) + beta_2 beta_1 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 3)), = 0) + beta^2_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 4)), = 0) = \
    & = beta_2 sigma^2 \
    rho_2 = frac(gamma_1, gamma_0) & = frac(beta_2, 1 + beta^2_1 + beta^2_2)
  $

+ Найдем третью автокорреляцию:
  $
    gamma_3 = Ex (y_t y_(t - 3)) & = Ex (epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2)) (epsilon_(t - 3) + beta_1 epsilon_(t - 4) + beta_2 epsilon_(t - 5)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 4)), = 0) + beta_2 underbrace(Ex (epsilon_t epsilon_(t - 5)), = 0) + \
    & wide + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 4)), = 0) + beta_1 beta_2 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 5)), = 0) + \
    & wide + beta_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 3)), = 0) + beta_2 beta_1 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 4)), = 0) + beta^2_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 4)), = 0) = \
    & = 0 \
    rho_3 = frac(gamma_1, gamma_0) & = 0
  $

Все последующие автокорреляции также равны 0.

Рассмотрим общий случай.
В числителе выражения для $phi_(s s)$
находится следующий определитель
$
  mat(
    delim: "|", rho_0, rho_1, rho_2, 0, 0, dots.h.c, 0, rho_(1);
    rho_1, rho_0, rho_1, rho_2, 0, dots.h.c, 0, rho_(2);
    rho_2, rho_1, rho_0, rho_1, rho_2, dots.h.c, 0, 0 ;
    0, rho_2, rho_1, rho_0, rho_1, dots.h.c, 0, 0 ;
    0, 0, rho_2, rho_1, rho_0, dots.h.c, 0, 0 ;
    dots.v, dots.v, dots.v, dots.v, dots.v, dots.down, dots.v, dots.v;
    0, 0, 0, 0, 0, dots.h.c, rho_0, 0 ;
    0, 0, 0, 0, 0, dots.h.c, rho_1, 0
  ) != 0
$
При этом, его значения с ростом $s$ стемятся к нулю.

== MA(q)

Проведя аналогичные расчеты для процесса MA(q) можно показать, что
$
  rho_s = cases(
    > 0\, & s = q,
    0\, & s > q
  )
$

При $s < q$ автокорреляции могут принимать любые допустимые значения.

Для PACF рассуджения аналогичны таковым для процесса MA(2).

== ARMA(1,1)

+ Найдем первую автокорреляцию:
  $
    Ex (y_t epsilon_t) & = Ex [ (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) epsilon_t ] = \
    & = a_1 underbrace(Ex y_(t - 1)epsilon_t, = 0) + underbrace(Ex epsilon^2_t, = sigma^2) + beta_1 underbrace(Ex epsilon_(t - 1)epsilon_t, = 0) = sigma^(2)
  $
  $
    Ex (y_t epsilon_(t - 1)) & = Ex [ (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) epsilon_(t - 1) ] = \
    & = a_1 underbrace(Ex (y_(t - 1)epsilon_(t - 1)), = sigma^2) + underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 1)), = sigma^2) = (a_1 + beta_1) sigma^(2)
  $
  $
    gamma_0 = Ex y^2_t & = Ex (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1))^2 = \
    & = a^2_1 underbrace(Ex y^2_(t - 1), = gamma_0) + a_1 underbrace(Ex (y_(t - 1)epsilon_t), = 0) + a_1 beta_1 underbrace(Ex (y_(t - 1)epsilon_(t - 1)), = sigma^2) + \
    & wide + a_1 underbrace(Ex (epsilon_t y_(t - 1)), = 0) + underbrace(Ex epsilon^2_t, = sigma^2) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + \
    & wide + a_1 beta_1 underbrace(Ex (epsilon_(t - 1)y_(t - 1)), = sigma^2) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_t), = 0) + beta^2_1 underbrace(Ex epsilon^2_(t - 1), = sigma^2) => \
    => gamma_0 & = frac(1 + beta^2_1 + 2 a_1 beta_1, 1 - a^2_1)sigma^(2)
  $
  $
    gamma_1 = Ex (y_t y_(t - 1)) & = Ex (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) (a_1 y_(t - 2) + epsilon_(t - 1) + beta_1 epsilon_(t - 2)) = \
    & = a^2_1 underbrace(Ex (y_(t - 1)y_(t - 2)), = gamma_1) + a_1 underbrace(Ex (y_(t - 1)epsilon_(t - 1)), = sigma^2) + a_1 beta_1 underbrace(Ex (y_(t - 1)epsilon_(t - 2)), = (a_1 + beta_1)sigma^2) + \
    & wide + a_1 underbrace(Ex (epsilon_t y_(t - 2)), = 0) + underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + \
    & wide + a_1 beta_1 underbrace(Ex (epsilon_(t - 1)y_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 1)), = sigma^2) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) = \
    & = a^2_1 gamma_1 + [ a_1 + a_1 beta_1 (a_1 + beta_1) + beta_1 ] sigma^2 => \
    => gamma_1 & = frac((1 + a_1 beta_1) (a_1 + beta_1), 1 - a^2_1) \
    rho_1 & = frac((1 + a_1 beta_1) (a_1 + beta_1), 1 + beta^2_1 + 2 a_1 beta_1)
  $

+ Найдем вторую автокорреляцию.
  $
    Ex (y_t epsilon_(t - 2)) & = Ex [ (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) epsilon_(t - 2) ] = \
    & = a_1 underbrace(Ex (y_(t - 1)epsilon_(t - 2)), = (a_1 + beta_1)sigma^2) + underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) = a_1 (a_1 + beta_1) sigma^(2)
  $
  $
    gamma_2 = Ex (y_t y_(t - 2)) & = Ex (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) (a_1 y_(t - 3) + epsilon_(t - 2) + beta_1 epsilon_(t - 3)) = \
    & = a^2_1 underbrace(Ex (y_(t - 1)y_(t - 3)), = gamma_2) + a_1 underbrace(Ex (y_(t - 1)epsilon_(t - 2)), = (a_1 + beta_1)sigma^2) + a_1 beta_1 underbrace(Ex (y_(t - 1)epsilon_(t - 3)), = a_1 (a_1 + beta_1)sigma^2) + \
    & wide + a_1 underbrace(Ex (epsilon_t y_(t - 3)), = 0) + underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + \
    & wide + a_1 beta_1 underbrace(Ex (epsilon_(t - 1)y_(t - 3)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) = \
    & = a^2_1 gamma_2 + [ a_1 (a_1 + beta_1) + a^2_1 beta_1 (a_1 + beta_1) ] sigma^2 => \
    => gamma_2 & = frac((1 + a_1 beta_1) (a_1 + beta_1), 1 - a^2_1) \
    rho_2 & = frac(a_1 (1 + a_1 beta_1) (a_1 + beta_1), 1 + beta^2_1 + 2 a_1 beta_1) = a_1 rho_(1)
  $

+ Найдем третью автокорреляцию.
  $
    Ex (y_t epsilon_(t - 3)) & = Ex [ (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) epsilon_(t - 3) ] = \
    & = a_1 underbrace(Ex (y_(t - 1)epsilon_(t - 3)), = a_1 (a_1 + beta_1)sigma^2) + underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) = a^2_1 (a_1 + beta_1) sigma^(2)
  $
  $
    gamma_3 = Ex (y_t y_(t - 3)) & = Ex (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) (a_1 y_(t - 4) + epsilon_(t - 3) + beta_1 epsilon_(t - 4)) = \
    & = a^2_1 underbrace(Ex (y_(t - 1)y_(t - 4)), = gamma_3) + a_1 underbrace(Ex (y_(t - 1)epsilon_(t - 3)), = a_1 (a_1 + beta_1)sigma^2) + a_1 beta_1 underbrace(Ex (y_(t - 1)epsilon_(t - 4)), = a^2_1 (a_1 + beta_1)sigma^2) + \
    & wide + a_1 underbrace(Ex (epsilon_t y_(t - 4)), = 0) + underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 4)), = 0) + \
    & wide + a_1 beta_1 underbrace(Ex (epsilon_(t - 1)y_(t - 4)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 4)), = 0) = \
    & = a^2_1 gamma_3 + [ a^2_1 (a_1 + beta_1) + a^3_1 beta_1 (a_1 + beta_1) ] sigma^2 => \
    => gamma_3 & = frac(a^2_1 (1 + a_1 beta_1) (a_1 + beta_1), 1 - a^2_1) \
    rho_3 & = frac(a^2_1 (1 + a_1 beta_1) (a_1 + beta_1), 1 + beta^2_1 + 2 a_1 beta_1) = a^2_1 rho_(1)
  $

Нетрудно показать, что
$ rho_s = a^(s - 1)_1 rho_1 $

Ручной расчет частных автокорреляций затруднен в силу одновременного присутствия AR и MA части.
Однако, наличие MA части приводит к тому, что значения PACF экспоненциально убывают.

== ARMA(p,q)

Процессы ARMA сочетают в себе "худшие" свойства автокорреляционной и частной автокорреляционной функций AR и MA процессов.
Для процесса $"ARMA"(p, q)$ наблюдается следующее:
- ACF начинает экспоненциально убывать, начиная с лага $q$; убывание может быть знакопеременным.
- PACF начинает экспоненциально убывать, начиная с лага $p$; убывание может быть знакопеременным.

Подытожим все вышесказанное в общей таблице.
#table(
  columns: (auto, auto, auto),
  stroke: none,
  align: (left, center, center),
  table.hline(),
  table.header([Процесс], [ACF], [PACF]),
  table.hline(),
  [WN], [$ rho_s = 0, forall s > 0 $], [$ phi_(s s) = 0, forall s > 0 $],
  [AR(1)],
  [$ rho_s = a^s_1, forall s > 0 $],
  [$
    phi_(s s) = cases(
      rho_(1)\, & s = 1,
      0\, & s > 1
    )
  $],
  [AR(2)],
  [$ rho_s = cases(frac(a_1, 1 - a_2)\, & s = 1, a_1 rho_(s - 1) + a_2 rho_(s - 2)\, & s > 1) $],
  [$
    phi_(s s) = cases(
      rho_(1)\, & s = 1,
      frac(rho_2 - rho^2_1, 1 - rho^2_1)\, & s = 2,
      0\, & s > 2
    )
  $],
  [AR(p)#footnote[В большинстве источников для $phi_(s s)$ принимают $s <= p$, но, строго говоря, требование неравенства нулю накладывается только на $phi_(p p)$.]],
  [$ rho_s = sum^p_(i = 1)a_i rho_(s - i) -> 0 $],
  [$
    phi_(s s) = cases(
      != 0\, & s = p,
      0\, & s > p
    )
  $],
  [MA(1)],
  [$
    rho_s = cases(
      frac(beta_1, 1 + beta^2_1)\, & s = 1,
      0\, & s > 1
    )
  $],
  [$
    phi_(s s) = cases(
      > 0\, & beta_(1)> 0 #text[и] s #text[нечетно],
      < 0\, & #text[иначе]
    ) -> 0
  $],
  [MA(2)],
  [$
    rho_s = cases(
      frac(beta_1 (1 + beta_2), 1 + beta^2_1 + beta^2_2)\, & s = 1,
      frac(beta_2, 1 + beta^2_1 + beta^2_2)\, & s = 2,
      0\, & s > 2
    )
  $],
  [$phi_(s s) -> 0$, возможно знакопеременно.],
  [MA(q)],
  [$
    rho_s = cases(
      != 0\, & s = q,
      0\, & s > q
    )
  $],
  [$phi_(s s) -> 0$, возможно знакопеременно.],
  [ARMA(p, q)],
  [$rho_s -> 0$, начиная с $q$], [$phi_(s s) -> 0$, начиная с $p$],
  table.hline(),
)

Обратите внимание на некоторую симметрию, присущую моделям AR и MA.

#bibliography("bibliography.bib")

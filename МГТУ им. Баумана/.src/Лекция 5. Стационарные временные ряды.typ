#import "@preview/theorion:0.6.0": *
#import "@preview/tschich:0.2.0": *
#import "@local/neat-document:0.1.0": *
#import "@local/trinkets:0.1.0": *
#set document(title: [Стационарные временные ряды])
#let margins = tschich-var(210mm, 297mm, 1 / 10)
#let font-size = 14pt
#show: neat-document.with(margin: margins, leading: 1em, fonts: (size: font-size))
#show: show-theorion
#show: trinkets

= Оператор запаздывания

/ Оператор запаздывания $L$: оператор, переводящий значения
временного ряда в момент времени $t$, в его значении в предыдущий
момент времени, то есть
$ L Y_t = Y_(t - 1) $

Оператор запаздывания обладает следующими свойствами:
+ $L c = c$;
+ _Дистрибутивный закон_: $(L^i + L^j) y_t = y_(t - i) + y_(t - j)$;
+ _Ассоциативный закон умножения_: $L^i L^j y_t = y_(t - i - j)$;
+ $L^(- s)y_t = y_(t + s)$;
+ Если $abs(a) < 1$, то $(1 + a L + a^2 L^2 + a^3 L^3 + ...) y_t = display(frac(1, 1 - a L))y_t$;
+ Если $abs(a) > 1$, то $(1 + (a L)^(- 1) + (a L)^(- 2) + (a L)^(- 3) + ...) y_t = - display(frac(a L, 1 - a L))y_t$.

= Стационарность

/ Сильная стационарность: независимость совместных распределений от сдвига во времени.
/ Слабая стационарность: выполнение следующих условий:
  - существование *математического ожидания*, его конечность и независимость от сдвига во времени;
  - существование *дисперсии*, ее конечность и независимость от сдвига во времени;
  - существование *автоковариаций*, их конечность и независимость от сдвига во времени; автоковариации должны зависеть *только* от взаимного расположения наблюдений.

*Автоковариационной функцией* случайного процесса $y_t$ называется функция $gamma (s) = gamma_s = Cov (y_t, y_(t - s))$.

Сильная стационарность не влечет слабую стационарность! Например, распределение Коши
$ f_X (x) = frac(1, pi gamma [ 1 + (frac(x - x_0, gamma))^2 ]) $ не имеет ни математического ожидания, ни дисперсии.

== Теорема Вольда

Любой слабостационарный недетерминированный случайный процесс $Y_t$ можно представить в виде
$ Y_t - mu = sum^oo_(i = 0) beta_i epsilon_(t - i) $
где $epsilon_t tilde.op "WN"(0, sigma^2_(epsilon))$.
При этом $sum^oo_(i = 0) abs(beta_i) < oo$.

== Слабая стационарность и обратимость MA процессов

*Процессом скользящего среднего* порядка $q$, или *процессом с ошибкой в форме скользящего среднего*, называется процесс вида
$ y_t = a_0 + epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2) + dots.h.c+ beta_q epsilon_(t - q) $
где $epsilon tilde.op "WN"(0,sigma^2)$.
Краткое обозначение:
$y_t tilde.op "MA"(q)$.

Случайный процесс называется *обратимым*, если его можно представить в виде процесса AR($oo$)
$ y_t = alpha_0 + sum^(oo)_(i = 1)alpha_i y_(t - i) + epsilon_t $
с абсолютно сходящимися весами $sum^(oo)_(i = 0)abs(alpha_i) < oo$.

Обычно определение обратимого процесса рассматривается применительно к процессам скользящего среднего.
Рассматривать его в отношении процесса авторегрессии *конечного* порядка нет особого смысла, так как процесс AR(p) является частным случаем процесса AR($oo$): достаточно принять $alpha_i = 0, forall i > p$.

=== MA(1)

Рассмотрим условия слабой стационарности процесса MA(1)
$ y_t = a_0 + epsilon_t + beta_1 epsilon_(t - 1) $

Найдем математическое ожидание процесса:
$
  Ex y_t & = Ex (a_0 + epsilon_t + beta_1 epsilon_(t - 1)) = \
         & = a_0 + underbrace(Ex epsilon_t, = 0) + beta_1 underbrace(Ex epsilon_(t - 1), = 0) = \
         & = a_0 = #text[const]
$

Найдем дисперсию процесса:
$
  Var y_t & = Ex (a_0 + epsilon_t + beta_1 epsilon_(t - 1) - a_0)^2 = \
  & = underbrace(Ex epsilon^2_t, = sigma^2) + 2 beta_1 underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta^2_1 underbrace(Ex epsilon^2_(t - 1), = sigma^2) = \
  & = (1 + beta^2_1) sigma^2 = #text[const]
$
Конечно, можно воспользоваться свойством дисперсии суммы независимых
случайных величин:
$
  Var y_t & = Var [ a_0 + epsilon_t + beta_1 epsilon_(t - 1) ] = \
          & = Var [ epsilon_t + beta_1 epsilon_(t - 1) ] = \
          & = Var epsilon_t + beta^2_1 Var epsilon_(t - 1) = \
          & = (1 + beta^2_1) sigma^(2)
$

Найдем автоковариации процесса.
$
  gamma_1 & = Cov (y_t,y_(t - 1)) = Ex [ (epsilon_t + beta_1 epsilon_(t - 1)) (epsilon_(t - 1) + beta_1 epsilon_(t - 2)) ] = \
  & = underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex epsilon^2_(t - 1), = sigma^2) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) = \
  & = beta_1 sigma^(2)
$
$
  gamma_2 & = Cov (y_t,y_(t - 2)) = Ex [ (epsilon_t + beta_1 epsilon_(t - 1)) (epsilon_(t - 2) + beta_1 epsilon_(t - 3)) ] = \
  & = underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) = \
  & = 0
$
Нетрудно показать, что $gamma_s = 0, forall s > 1$.

Для получения представления в форме AR($oo$) нужно к обеим частям процесса применить оператор, обратный к $(1 + beta_1 L)$:
$ (1 + beta_1 L)^(- 1)y_t = (1 + beta_1 L)^(- 1)a_0 + epsilon_t $

Для того, чтобы такой оператор существовал, необходимо выполнение условия $|beta_1| < 1$; тогда, воспользовавшись пятым свойством лагового оператора, получим
$
  (1 + beta_1 L)^(- 1)y_t & = (1 - beta_1 L + beta^2_1 L^2 - beta^3_1 L^3 + beta^4_1 L^(4)...) y_t = \
                          & = y_t - beta_1 y_(t - 1) + beta^2_1 y_(t - 2) + dots.h.c= \
                          & = y_t + sum^(oo)_(i = 1)(- 1)^i beta^i_1 y_(t - i) \
  (1 + beta_1 L)^(- 1)a_0 & = (1 - beta_1 L + beta^2_1 L^2 + dots.h.c) a_0 = \
                          & = a_0 - beta_1 a_0 + beta^2_1 a_0 + dots.h.c= \
                          & = a_0 (1 + sum^(oo)_(i = 1)(- 1)^i beta^i_1) = \
                          & = frac(a_0, 1 + beta_1)
$

Отсюда искомое представление имеет вид
$ y_t = frac(a_0, 1 + beta_1) - sum^(oo)_(i = 1)beta^i_1 y_(t - i) + epsilon_t $

=== MA(q)

Рассмотрим условия слабой стационарности процесса MA(q)
$ y_t = a_0 + epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2) + dots.h.c+ beta_q epsilon_(t - q) $

Математическое ожидание процесса равно:
$
  Ex y_t = Ex (a_0 + epsilon_t + sum^q_(i = 1)beta_i epsilon_(t - i)) = a_0 + Ex epsilon_t + sum^q_(i = 1)beta_i Ex epsilon_(t - i) = a_0 = #text[const]
$

Дисперсия процесса равна:
$
  Var y_t & = Var [ a_0 + epsilon_t + sum^q_(i = 1)beta_i epsilon_(t - i) ] = \
          & = Var epsilon_t + sum^q_(i = 1)beta^2_q Var epsilon_(t - i) = (1 + sum^q_(i = 1)beta^2_i) sigma^(2)
$

Авторегрессии процесса равны:
$
  gamma_s = cases(
    (sum^(q - s)_(i = 0)beta_i beta_(i + s)) sigma^(2)\, & s = overline(1\,q - 1),
    beta_q sigma^(2)\, & s = q,
    0\, & s > q
  )
$

Для получения представления в форме AR($oo$) нужно к обеим частям процесса применить оператор, обратный к $(1 + beta_1 L + beta_2 L^2 + ...+ beta_q L^q)$.
Отметим, что лаговый полином порядка $q$, который нам следует обратить, можно разложить на полиномы первого порядка следующим образом:
$
  1 + beta_1 L + beta_2 L^2 + ...+ beta_q L^q = (1 - pi_1 L) (1 - pi_2 L) ...
  (1 - pi_q L)
$
где $pi_1,..., pi_q in CC$.
Такое разложение всегда существует согласно Основной теореме алгебры; более того, оно единственно.

С учетом этого обращение полинома сводится к обращению $q$ полиномов первой степени.
Условие существования обратных полиномов такое же:
$abs(pi_i) < 1, i = overline(#[$1, q$])$.

Выражение для свободного члена получается аналогично процессу MA(1) путем последовательного применения к $a_0$ всех обратных операторов.

=== MA(inf)

Рассмотрим условия слабой стационарности процесса $"MA"(oo)$
$ y_t = a_0 + epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2) + ... $

Математическое ожидание процесса равно:
$ Ex y_t = a_0 $

Дисперсия процесса равна:
$ Var y_t = sigma^2 sum^(oo)_(i = 0)beta^2_i $

Авторегрессии процесса равны:
$ gamma_s = sigma^2 sum^(oo)_(i = 0)beta_i beta_(i + s) $

Как можно понять из предыдущих результатов, процесс скользящего среднего *конечного* порядка *всегда* слабо стационарен.
Для бесконечного MA процесса это не так: для существования вторых моментов процесса необходимо выполнение следующего условия
$ sum^(oo)_(i = 0)abs(beta_i) < oo $

== Слабая стационарность AR процессов

*Процессом авторегрессии* порядка $p$ называется процесс вида
$ y_t = a_0 + epsilon_t + a_1 y_(t - 1) + a_2 y_(t - 2) + dots.h.c+ a_p y_(t - p) + epsilon_t $
где $epsilon tilde.op "WN"(0,sigma^2)$.
Краткое обозначение:
$y_t tilde.op "AR"(p)$.

Рассмотрим модель ARMA(p,q)
$
  y_t = a_0 + a_1 y_(t - 1) + a_2 y_(t - 2) + dots.h.c+ a_p y_(t - p) + epsilon_t + beta_1 epsilon_(t - 1) + dots.h.c+ beta_q epsilon_(t - q)
$
Если перенести все запаздывающие значения $y_t$ в левую часть, то получим, что данная модель представляет собой неоднородное разностное
уравнение порядка $p$
$
  underbrace(y_t - a_1 y_(t - 1) - a_2 y_(t - 2) + dots.h.c- a_p y_(t - p), #text[однородное разностное уравнение]) = underbrace(a_0 + epsilon_t + beta_1 epsilon_(t - 1) + dots.h.c+ beta_q epsilon_(t - q), #text[неоднородная] #text[часть])
$

Если решить это уравнение, что можно найти выражение для $y_t$ не в рекурсивной форме.
Напомним, что из курса дифференциальных и разностных уравнений (см., например, #cite(<romanko2012kurs>)) следует, что общее решение такого уравнения равно сумме решений однородной части и какого-либо частного решения
$ y_t = y^(#text[однородное])_t + y^(#text[частное])_t $

Подробно рассмотрим решение этого уравнения на примере авторегрессии первого порядка.
Для авторегрессий большего порядка запишем только общее решение.

=== AR(1)

Рассмотрим модель
$ y_t = a_0 + a_1 y_(t - 1) + epsilon_t $ <eq-4.1>
Найдем решение однородного уравнения
$ y_t = a_1 y_(t - 1) $

Очевидно, что тривиальное решение $y_t equiv 0$ нас не сильно устраивает, хотя оно, несомненно, является решением данного однородного уравнения.
Теория говорит нам, что имеет смысл искать нетривиальное решение данного уравнения в виде $y_t = lambda^t$; отметим, что таким образом мы ищем общее решение с точностью до ненулевого множителя.

Подставим это выражение в уравнение
$ lambda^t = a_1 lambda^(t - 1) $
Поскольку нас интересует нетривиальное решение системы, мы можем смело предполагать $lambda != 0$ и не менее смело сократить данное выражение на $lambda^(t - 1)$ и получить
#eqnum[$ lambda = a_1, $] <eq-4.2>
то есть решение однородной части имеет вид
$ y^(#text[однородное])_t = C a^t_1 $

Для того, чтобы это решение было устойчивым, необходимо, чтобы выполнялось условие $|a_(1)| < 1$, что совпадает с условием существования представления процесса AR(1) в виде MA($oo$).

Найдем теперь частное решение неоднородного уравнения.
Воспользуемся *методом неопределенных коэффициентов*; для этого будем искать частное решение в виде
$ y^(#text[частное])_t = alpha + sum^(oo)_(i = 0)beta_i epsilon_(t - i) $

Подставив это выражение в уравнение @eq-4.1, получим
$
  alpha + sum^(oo)_(i = 0)beta_i epsilon_(t - i) = a_0 + a_1 alpha + a_1 sum^(oo)_(i = 0)beta_i epsilon_(t - i - 1) + epsilon_(t)
$
Приравняв коэффициенты при подобных слагаемых, получим
$
   alpha & = a_0 + a_1 alpha, => &  alpha & = frac(a_0, 1 - a_1) \
  beta_0 & = 1 \
  beta_1 & = a_1 beta_0, =>      & beta_1 & = a_1 \
  beta_2 & = a_1 beta_1, =>      & beta_2 & = a^2_1 \
     ...
$

Итого
$ y^(#text[частное])_t = frac(a_0, 1 - a_1) + sum^(oo)_(i = 0)a^i_1 epsilon_(t - i), $
а полное решение
$ y_t = C_0 a^t_1 + frac(a_0, 1 - a_1) + sum^(oo)_(i = 0)a^i_1 epsilon_(t - i) $

То же самое решение можно получить методом *вариации константы*.
Будем искать общее решение в виде
#eqnum[$ y_t = C_t a^t_1, $] <eq-4.3>
то есть будем рассматривать постоянный множитель как функцию от времени.
Подставив это выражение в уравнение @eq-4.1, получим
$
  C_t a^t_1 & = a_0 + a_1 C_(t - 1)a^(t - 1)_1 + epsilon_t = a_0 + C_(t - 1)a^t_1 + epsilon_t \
        C_t & = frac(a_0, a^t_1) + C_(t - 1) + frac(epsilon_t, a^t_1)
$
Начнем рекурсивно подставлять $C_t$ в эту формулу; получим
$
  C_t & = frac(a_0, a^t_1) + C_(t - 1) + frac(epsilon_t, a^t_1) = \
  & = frac(a_0, a^t_1) + frac(a_0, a^(t - 1)_1) + C_(t - 2) + frac(epsilon_(t - 1), a^(t - 1)_1) + frac(epsilon_t, a^t_1) = \
  &...\
  & = sum^(t - 1)_(i = 0)frac(a_0, a^(t - i)_1) + C_0 + sum^(t - 1)_(i = 0)frac(epsilon_(t - i), a^(t - i)_1)
$
Подставим полученное выражение в @eq-4.2
$
  y_t & = (sum^(t - 1)_(i = 0)frac(a_0, a^(t - i)_1) + C_0 + sum^(t - 1)_(i = 0)frac(epsilon_(t - i), a^(t - i)_1)) a^t_1 = \
  & = underbrace(C_0 a^t_1, ->_(t -> oo)0) + underbrace(sum^(t - 1)_(i = 0)a_0 a^i_1, ->_(t -> oo)frac(a_0, 1 - a_1)) + sum^(t - 1)_(i = 0)a^i_1 epsilon_(t - i)
$

Несмотря на отличия в пределах суммирования, общий вид решения совпадает.
Легко увидеть, что $C_0 = y_0$.

Видно, что если $abs(a_1) < 1$, то
- $y_0 a^t_i limits(->)_(t -> oo) 0$;
- свободный член стремится к значению математического ожидания процесса;
- последнее слагаемое — процесс MA($oo$) — является слабостационарным процессом.

То есть общее решение разностного уравнения стремится к представлению Вольда процесса AR(1).

Отметим, что выражение @eq-4.3 по своей сути является результатом подстановки $lambda^t$ вместо $y_t$ в однородную часть разностного уравнения.

Данное уравнение называется *характеристическим уравнением*.
Условием стабильности решения разностного уравнения является выполнение следующего условия
$ abs(lambda_i) < 1, forall i $
где $lambda_i$ — корни характеристического уравнения.

В рассматриваемом примере процесса AR(1) такой корень один, и он равен $a_1$, откуда и следует рассмотренное ранее условие $abs(a_1) < 1$.
Если же мы подставили бы вместо $L$ какую-либо переменную $z$ и решали бы уравнение $A(z) = 0$, то мы получили бы *обратное характеристическое уравнение*.
Условием стабильности для обратного характеристического уравнения является
$ abs(z_i) > 1, forall i $
где $z_i$ — корни обратного характеристического уравнения.

Корни этих двух уравнений связаны следующим образом: для любого корня $lambda_i$ характеристического уравнения существует соответствующий ему корень обратного характеристического уравнения $z_i$ такой, что $lambda_i = frac(1, z_i)$.

Наконец, отметим, что стабильность решения разностного уравнения равносильна слабой стационарности соответствующего авторегрессионного процесса.

=== AR(p)

Рассмотрим модель
$ y_t = a_0 + a_1 y_(t - 1) + a_2 y_(t - 2) + dots.h.c+ a_p y_(t - p) + epsilon_t $
Соответствующее характеристическое уравнение представляет собой полином степени $p$, который имеет, соответственно, $p$ корней $lambda_i in CC$.
Условие стационарности не меняется: $abs(lambda_i) < 1, forall i = overline(#[$1, p$])$.

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
- $r_s limits(->)_(T -> oo) rho_s$.

Автокорреляция, являясь важной характеристикой случайного процесса, обладает определенным недостатком.
Давайте рассмотрим процесс AR(2)
$ y_t = a_1 y_(t - 1) + a_2 y_(t - 2) + epsilon_t $

Если мы найдем значение $rho_2$, то оно, строго говоря, будет показывать нечто среднее между "прямой" корреляцией $y_t$ и $y_(t - 2)$ и "опосредованной" корреляцией, возникающей из-за того, что $y_t$ зависит от величины $y_(t - 1)$, которая, в свою очередь, зависит от $y_(t - 2)$.
#align(center)[#image("fig/PACF.pdf", width: 60%)]

Для того, чтобы обойти этот эффект, используют *частные автокорреляции*.
*Частная автокорреляция* порядка $s$ определяется как значение коэффициента $phi_(s s)$ в следующей модели AR(s):
#eqnum[$ y_t = phi_(0 s) + phi_(s 1)y_(t - 1) + phi_(s 2)y_(t - 2) + dots.h.c+ phi_(s s)y_(t - s) + nu_t $] <eq-6.1>

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
  macron(y)_t & = y_t - mu = a_0 - frac(a_0, 1 - sum^p_(i = 1)a_i) + a_1 y_(t - 1) + dots.h.c+ a_p y_(t - p) + epsilon_t = \
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
    gamma_0 & = Ex y^2_t = Ex (epsilon_t + beta_1 epsilon_(t - 1))^2 = \
    & = underbrace(Ex epsilon^2_t, = sigma^2) + 2 beta_1 underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta^2_1 underbrace(Ex epsilon^2_(t - 1), = sigma^2) = \
    & = (1 + beta^2_1)sigma^2 \
    gamma_1 & = Ex (y_t y_(t - 1)) = Ex (epsilon_t + beta_1 epsilon_(t - 1)) (epsilon_(t - 1) + beta_1 epsilon_(t - 2)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 1)), = sigma^2) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) = \
    & = beta_1 sigma^2 \
    rho_1 & = frac(gamma_1, gamma_0) = frac(beta_1, 1 + beta^2_1)
  $

  Обратим внимание, что для процесса MA(1) первая автокорреляция не превышает по абсолютному значению $0.5$!

+ Найдем вторую автокорреляцию:
  $
    gamma_2 & = Ex (y_t y_(t - 2)) = Ex (epsilon_t + beta_1 epsilon_(t - 1)) (epsilon_(t - 2) + beta_1 epsilon_(t - 3)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) = \
    & = 0 \
    rho_2 & = frac(gamma_2, gamma_0) = 0
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
    gamma_0 & = Ex y^2_t = Ex (epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2))^2 = \
    & = underbrace(Ex epsilon^2_t, sigma^2) + beta^2_1 underbrace(Ex epsilon^2_(t - 1), sigma^2) + beta^2_2 underbrace(Ex epsilon^2_(t - 2), sigma^2) + \
    & wide + 2 beta_1 underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + 2 beta_2 underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + 2 beta_1 beta_2 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) = \
    & = (1 + beta^2_1 + beta^2_2) sigma^(2)
  $
  $
    gamma_1 & = Ex (y_t y_(t - 1)) = Ex (epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2)) (epsilon_(t - 1) + beta_1 epsilon_(t - 2) + beta_2 epsilon_(t - 3)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_2 underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + \
    & wide + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 1)), = sigma^2) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) + beta_1 beta_2 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) + \
    & wide + beta_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 1)), = 0) + beta_2 beta_1 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 2)), = sigma^2) + beta^2_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 3)), = 0) = \
    & = (beta_1 + beta_1 beta_2) sigma^2 \
    rho_1 & = frac(gamma_1, gamma_0) = frac(beta_1 + beta_1 beta_2, 1 + beta^2_1 + beta^2_2)
  $

+ Найдем вторую автокорреляцию:
  $
    gamma_2 & = Ex (y_t y_(t - 2)) = Ex (epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2)) (epsilon_(t - 2) + beta_1 epsilon_(t - 3) + beta_2 epsilon_(t - 4)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 2)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + beta_2 underbrace(Ex (epsilon_t epsilon_(t - 4)), = 0) + \
    & wide + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 2)), = 0) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) + beta_1 beta_2 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 4)), = 0) + \
    & wide + beta_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 2)), = sigma^2) + beta_2 beta_1 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 3)), = 0) + beta^2_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 4)), = 0) = \
    & = beta_2 sigma^2 \
    rho_2 & = frac(gamma_1, gamma_0) = frac(beta_2, 1 + beta^2_1 + beta^2_2)
  $

+ Найдем третью автокорреляцию:
  $
    gamma_3 & = Ex (y_t y_(t - 3)) = Ex (epsilon_t + beta_1 epsilon_(t - 1) + beta_2 epsilon_(t - 2)) (epsilon_(t - 3) + beta_1 epsilon_(t - 4) + beta_2 epsilon_(t - 5)) = \
    & = underbrace(Ex (epsilon_t epsilon_(t - 3)), = 0) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 4)), = 0) + beta_2 underbrace(Ex (epsilon_t epsilon_(t - 5)), = 0) + \
    & wide + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 3)), = 0) + beta^2_1 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 4)), = 0) + beta_1 beta_2 underbrace(Ex (epsilon_(t - 1)epsilon_(t - 5)), = 0) + \
    & wide + beta_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 3)), = 0) + beta_2 beta_1 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 4)), = 0) + beta^2_2 underbrace(Ex (epsilon_(t - 2)epsilon_(t - 4)), = 0) = \
    & = 0 \
    rho_3 & = frac(gamma_1, gamma_0) = 0
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
    gamma_0 & = Ex y^2_t = Ex (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1))^2 = \
    & = a^2_1 underbrace(Ex y^2_(t - 1), = gamma_0) + a_1 underbrace(Ex (y_(t - 1)epsilon_t), = 0) + a_1 beta_1 underbrace(Ex (y_(t - 1)epsilon_(t - 1)), = sigma^2) + \
    & wide + a_1 underbrace(Ex (epsilon_t y_(t - 1)), = 0) + underbrace(Ex epsilon^2_t, = sigma^2) + beta_1 underbrace(Ex (epsilon_t epsilon_(t - 1)), = 0) + \
    & wide + a_1 beta_1 underbrace(Ex (epsilon_(t - 1)y_(t - 1)), = sigma^2) + beta_1 underbrace(Ex (epsilon_(t - 1)epsilon_t), = 0) + beta^2_1 underbrace(Ex epsilon^2_(t - 1), = sigma^2) => \
    => gamma_0 & = frac(1 + beta^2_1 + 2 a_1 beta_1, 1 - a^2_1)sigma^(2)
  $
  $
    gamma_1 & = Ex (y_t y_(t - 1)) = Ex (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) (a_1 y_(t - 2) + epsilon_(t - 1) + beta_1 epsilon_(t - 2)) = \
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
    gamma_2 & = Ex (y_t y_(t - 2)) = Ex (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) (a_1 y_(t - 3) + epsilon_(t - 2) + beta_1 epsilon_(t - 3)) = \
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
    gamma_3 & = Ex (y_t y_(t - 3)) & = Ex (a_1 y_(t - 1) + epsilon_t + beta_1 epsilon_(t - 1)) (a_1 y_(t - 4) + epsilon_(t - 3) + beta_1 epsilon_(t - 4)) = \
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
      > 0\, & beta_(1) > 0 #text[и], & s #text[нечетно],
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

= Процедура Бокса-Дженкинса

+ Определение порядка интегрированности процесса.
+ Определение порядков ARMA процесса.
+ Оценивание модели процесса.
+ Диагностика оцененной модели.

== Четвертый этап

=== Автокорреляция в остатках

Для тестирования гипотезы об отсутствии автокорреляции в остатках модели можно использовать тесты Бокса-Пирса и Льюнг-Бокса, но важно помнить, что у тестовых статистик меняется распределение: при применении тестов к остаткам оценивания модели ARMA(p,q) статистики будут иметь асимптотическое распределение $chi^2 (s - p - q)$.

*Тест Бройша-Годфри* представляет собой тест множителей Лагранжа.
Процедура тестирования включает в себя два шага:

+ При помощи МНК оцениваем регрессию
  $ y_t = beta_0 + beta_1 x_(1 t) + dots.h.c+ beta_k x_(k t) + epsilon_t, $
  где $x_(i t)$ — какие-то интересующие нас регрессоры, в том числе лаги $y_t$.
  После оценки регрессии получаем остатки $hat(epsilon)_t$.
+ Оцениваем регрессию
  $
    hat(epsilon)_t = alpha_0 + alpha_1 x_(1 t) + dots.h.c+ alpha_k x_(k t) + gamma_1 hat(epsilon)_(t - 1) + dots.h.c+ gamma_s hat(epsilon)_(t - s) + nu_(t)
  $
  Запоминаем значение *нескорректированного (!)* $R^2$.
  Для тестирования набора гипотез
  $
    H_0 & : gamma_1 = gamma_2 = dots.h.c= gamma_s = 0 \
    H_a & : sum^s_(i = 1)gamma^2_(i)> 0
  $
  воспользуемся тестовой статистикой
  $ T R^2 tilde.op chi^2 (s) $
  Если остатки имеют нормальное распределение, то можно проводить тестирование с помощью F-статистики.

=== Нормальность остатков

Для этого можно воспользоваться *тестом Харке-Бера*, в котором рассматриваются следующие гипотезы
$
  H_0 & : epsilon tilde.op N (0,sigma^2) \
  H_a & : epsilon tilde.not N (0,sigma^2)
$
Строго говоря в тесте проверяется не нормальность, а совпадение коэффициентов асимметрии и эксцесса у распределения ошибок со значениями, присущими нормальному распределению (0 и 3 соответственно).

Статистика имеет вид
$
  J B & = T/6[ hat(S) + 1/4(hat(K) - 3)^2 ] limits(tilde)^(a s y .) chi^2 (2) \
  hat(S) & = frac(1/T sum^T_(t = 1)(hat(epsilon)_t - hat(mu)_(hat(epsilon)))^3, hat(sigma)^3_(hat(epsilon)))quad #text[выборочный коэффициент асимметрии] \
  hat(K) & = frac(1/T sum^T_(t = 1)(hat(epsilon)_t - hat(mu)_(hat(epsilon)))^4, hat(sigma)^4_(hat(epsilon)))quad #text[выборочный коэффициент эксцесса]
$

Тест Харке-Бера очень чувствителен к наличию выбросов в данных.

=== Условная гетероскедастичность

*ARCH-тест*, так же как и тест Бройша-Годфри, представляет собой LM-тест.
Процедура тестирования включает в себя два шага:

+ При помощи МНК оцениваем регрессию
  $ y_t = beta_0 + beta_1 x_(1 t) + dots.h.c+ beta_k x_(k t) + epsilon_t, $
  где $x_(i t)$ — какие-то интересующие нас регрессоры, в том числе лаги $y_t$.
  После оценки регрессии получаем остатки $hat(epsilon)_t$.
+ Оцениваем регрессию
  $
    hat(epsilon)^2_t = alpha_0 + gamma_1 hat(epsilon)^2_(t - 1) + dots.h.c+ gamma_s hat(epsilon)^2_(t - s) + nu_(t)
  $
  Запоминаем значение *нескорректированного (!)* $R^2$.
+ Для тестирования набора гипотез
  $
    H_0 & : gamma_1 = gamma_2 = dots.h.c= gamma_s = 0 \
    H_a & : sum^s_(i = 1)gamma^2_(i)> 0
  $
  воспользуемся тестовой статистикой
  $ T R^2 tilde chi^2 (s) $

Если остатки имеют нормальное распределение, то можно проводить тестирование с помощью F-статистики.
Если нулевая гипотеза не отвергается, то ошибки модели условно гомоскедастичны; в противном случае они условно гетероскедастичны.

#bibliography("bibliography.bib")

#show: appendix

= Корни характеристического уравнения и разложение лагового полинома AR(p)

Перепишем процесс AR(p) через лаговый оператор
$ underbrace([ 1 - a_1 L - a_2 L^2 - dots.h.c - a_p L^p ], A(L)) y_t = a_0 + epsilon_t $

Покажем, что $lambda_i = pi_i$ для AR(p).
Здесь $pi_i$ --- коэффициенты разложения лагового полинома
$
  A(L) = [1 - pi_1 L][1 - pi_2 L] ... [1 - pi_p L]
$

Заменив $L$ на $z$ и приравняв $A(z) = 0$, получим _обратное_ характеристическое уравнение
$ 1 - a_1 z - a_2 z^2 - dots.h.c - a_p z^p = 0, #<eq-p4.1.1> $
которое, согласно основной теореме алгебры, имеет следующее _единственное_ разложение
$ A(z) = -a_p (z - z_1) (z - z_2) dots.h.c (z - z_p) = 0 = bigast $

Нетрудно показать, что $1 = -a_p (-1)^p product_(i = 1)^p z_i$ и $-a_p = (-1)^p product_(i = 1)^p frac(1, z_i)$.
Отсюда
$
  bigast & = (-1)^p (product^p_(i = 1) frac(1, z_i)) (z - z_1) (z - z_2) dots.h.c (z - z_p) = nonumber \
         & = (-1)^p product^p_(i = 1) frac(z - z_i, z_i) = nonumber \
         & = product^p_(i = 1) autolr(1 - underbrace(frac(1, z_i), pi_i) z, delim: #("[", "]")) #<eq-p4.1.2>
$

Заменой $z = 1/lambda$ выражение @eq-p4.1.2 представимо в виде
$ A(z) = frac(1, lambda^p) product^p_(i = 1)(lambda - pi_i) #<eq-p4.1.3> $

Той же подстановкой из выражения @eq-p4.1.1 получим
$ A(z) = frac(1, lambda^p)(lambda^p- a_1 lambda^(p-1) - dots.h.c - a_p) = bigast $
где выражение в скобках представляет собой _прямое_ характеристическое уравнение, имеющее разложение
$ bigast = frac(1, lambda^p) product^p_(i = 1) (lambda - lambda_i) #<eq-p4.1.4> $

Выражения @eq-p4.1.3 и @eq-p4.1.4 совпадают с точностью до обозначений.
В силу того, что @eq-p4.1.4 получается из $A(z)$ заменой переменных, то обратной заменой получим @eq-p4.1.2 с точностью до множителей $lambda_i,pi_i$; поскольку разложение единственно, то $lambda_i = pi_i$.

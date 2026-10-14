#import "@preview/theorion:0.6.0": *
#import "@preview/tschich:0.2.0": *
#import "@local/neat-document:0.1.0": *
#import "@local/trinkets:0.1.0": *
#set document(title: [Модели дискретного выбора])
#show: neat-document.with(margin: tschich-var(210mm, 297mm, 1 / 12))
#show: show-theorion

Рассмотрим модели, в которых объясняемая переменная может принимать два значения (или несколько значений).

= Линейная модель вероятности

Простейшая модель, имеющая вид
$ vv(y) = vv(X) vv(beta) + vv(u) $
где $vv(y)$ принимает значения $0$ или $1$, и $Ex (vv(u)) = vv(0)$.

Условное математическое ожидание $i$-го наблюдения равно
$ Ex (vv(y)_i mid(bar.v) vv(x)_i) = vv(x)^prime_i vv(beta) $
Данное математическое ожидание следует интерпретировать как вероятность $Pr (y_i = 1)$ #cite(<pierse2021limited>).

В этом месте возникает первая проблема линейной модели:
#box(width: 100%, inset: 0.5em, stroke: 0.2em)[
  Линейная модель вероятности не может гарантировать того, что условное математическое ожидание будет ограничено между $0$ и $1$!
]

Более того, так как $y_i$ принимает только два значения: $0$ и $1$ — то ошибка $u_i$ также принимает только два значения
$
  u_i = cases(
    -vv(x)^prime_i vv(beta) & " если " y_i = 0,
    1 - vv(x)^prime_i vv(beta) & " если " y_i = 1
  )
$

Из условия $Ex (u_i) = 0$ следует, что
$
  Pr (- vv(x)^prime_i vv(beta)) + (1 - p) (1 - vv(x)^prime_i vv(beta)) = 0 \
  1 - p - vv(x)^prime_i vv(beta) = 0 \ Pr= 1 - vv(x)^prime_i vv(beta)
$
то есть $Pr (- vv(x)_i^prime vv(beta)) = 1 - vv(x)_i^prime vv(beta)$ и $Pr (1 - vv(x)_i^prime vv(beta)) = vv(x)_i^prime vv(beta)$.
Отсюда легко увидеть, что
$
  Var[u_i] & = sum u^2_i Pr (u_i) = \
  & = (- vv(x)^prime_i vv(beta))^2 (1 - vv(x)^prime_i vv(beta)) + (1 - vv(x)^prime_i vv(beta))^2 vv(x)^prime_i vv(beta) = \
  & = vv(x)^prime_i vv(beta) (1 - vv(x)^prime_i vv(beta))
$
то есть дисперсия ошибки _гетероскедастична_, и МНК оценки будут неэффективными.

Если рассмотреть подобранные значения
$ hat(y)_i = vv(x)^prime_i hat(vv(beta)) $
то оценкой дисперсии $Var[u_i]$ будет $hat(sigma)^2_i = hat(y)_i (1 - hat(y)_i)$.
Используя эти оценки, можно построить матрицу
$
  vv(L) = mat(
    delim: "(", frac(1, hat(sigma)_1);
    , dots.down;
    , , frac(1, hat(sigma)_N)
  )
$
которую следует использовать для доступного ВМНК.
#box(width: 100%, inset: 0.5em, stroke: 0.2em)[
  Как уже было сказано ранее, линейная модель вероятности не может гарантировать того, что условное математическое ожидание будет ограничено между $0$ и $1$.
  Нетрудно убедиться, что в этом случае $hat(sigma)^2_i< 0$.
]

= Пробит, логит и другие подобные модели

Давайте рассмотрим немного другой подход к проблеме.
Предположим, что у нас есть некоторая _латентная_#footnote[Скрытая и ненаблюдаемая.] переменная $vv(y)^(*)$, для которой построена модель регрессии
$ vv(y)^(*) = vv(X) vv(beta) + vv(u) $ <eq-ch5-y-star>
При этом _наблюдаемая_ переменная $vv(y)$ связана с ней следующим простым соотношением
#eqnum[$
  vv(y) = cases(
    1 " если " vv(y)^(*) > 0\,,
    0 " если " vv(y)^(*)lt.eq.slant 0
  ) = bb(I) [ vv(y)^(*)> 0 ]
$] <eq-ch5-y>

Как отмечает Пирс #cite(<pierse2021limited>), у нас нет никаких оснований считать, что $vv(u)$ и, соответственно, $vv(y)^(*)$ не имеют какого-либо непрерывного распределения $F$.
Более того, так как формула @eq-ch5-y не зависит от масштаба переменной $vv(y)^(*)$, то мы можем принять любую удобную для нас нормировку, например $Var[u_i] = 1$.

Из формул @eq-ch5-y-star и @eq-ch5-y следует, что
$
  P_i = Pr (y_i = 1) & = Pr (y^(*)_i> 0) = \
                     & = Pr (u_i> - vv(x)^prime_i vv(beta)) = \
                     & = 1 - Pr (u_i lt.eq.slant - vv(x)^prime_i vv(beta)) = \
                     & = 1 - F(- vv(x)^prime_i vv(beta))
$
где $F$ — кумулятивная функция распределения#footnote[Мы ее знаем как просто функцию распределения.]
$ F (x) = integral^x_(- infinity)f(u) dif u $

Если в качестве распределения ошибок используется какое-либо симметричное распределение — нормальное, Стьюдента и так далее — то последняя формула упрощается до
$ Pr (y_i = 1) = F(vv(x)^prime_i vv(beta)) $

В свете сказанного выше понятно, что $y_i$ имеет биномиальное распределение, а функция правдоподобия будет иметь вид
$ L (vv(y)) = product_(y_i = 1)P_i product_(y_i = 0)(1 - P_i) $

Если $vv(u)$ будет иметь #cite(<nosko2021ekonometrikaa>)
- нормальное распределение с
  $ Phi (x) = frac(1, sqrt(2 pi))integral^x_(- infinity)exp { - frac(u^2, 2) } dif u, $
  то полученная модель будет называться _пробит моделью_;
- логистическое распределение с
  $ Lambda (x) = frac(exp { x }, 1 + exp { x }) = frac(1, 1 + exp { - x }), $
  то полученная модель будет называться _логит моделью_;
- распределение Гомпертца c
  $ G (x) = 1 - exp { - e^x }, $
  то полученная модель будет называться _гомпит моделью_.

== Интерпретация коэффициентов

Так как описанный выше модели являются нелинейными по параметрам, то их интерпретация в традиционном смысле будет затруднительной.

Рассмотрим пробит модель
$ P_i = Pr (y_i = 1) = Phi (vv(x)^prime_i vv(beta)) $

Посчитаем предельный эффект от изменения регрессора $x_(i j)$
$ frac(partial P_i, partial x_(i j)) = phi(vv(x)^prime_i vv(beta)) beta_j $

Нетрудно увидеть, что, во-первых, он не равен $beta_j$ как в случае с линейными моделями, и, во-вторых, предельный эффект зависит от того, в какой точке он считается.

Рассмотрим теперь логит модель
$ P_i = Pr (y_i = 1) = Lambda (vv(x)^prime_i vv(beta)) $
которая несложными преобразованиями приводится к следующему виду
$ ln frac(P_i, 1 - P_i) = vv(x)^prime_i vv(beta) $
то есть натуральный логарифм _соотношения шансов_ в для логит модели линеен относительно $vv(x)_i$.
Предельный эффект же для логит модели равен
$
  frac(partial P_i, partial x_(i j)) & = lambda (vv(x)^prime_i vv(beta)) beta_j = frac(partial, partial x_(i j))[ frac(exp { vv(x)^prime_i vv(beta) }, 1 + exp { vv(x)^prime_i vv(beta) }) ] = \
  & = frac(exp { vv(x)^prime_i vv(beta) } (1 + exp { vv(x)^prime_i vv(beta) }) beta_j - exp { vv(x)^prime_i vv(beta) }^2 beta_j, (1 + exp { vv(x)^prime_i vv(beta) })^2) = \
  & = frac(exp { vv(x)^prime_i vv(beta) }, (1 + exp { vv(x)^prime_i vv(beta) })^2)beta_j = frac(1, 1 + exp { vv(x)^prime_i vv(beta) })(1 - frac(1, 1 + exp { vv(x)^prime_i vv(beta) })) beta_j = \
  & = P_i (1 - P_i) beta_j
$

== Меры качества моделей

Для рассмотрения мер качества подгонки введем, для начала, способ расчета подобранных значений.
Носко отмечает #cite(<nosko2021ekonometrikaa>), что наиболее естественным является
$
  hat(y_i) = cases(
    1\, & F(vv(x)^prime_i vv(beta)) > 1/2\,,
    0\, & F(vv(x)^prime_i vv(beta)) lt.eq.slant 1/2
  )
$

Рассмотрим классическую формулу для коэффициента детерминации
$ R^2 = 1 - frac(sum^N_(i = 1)(y_i - hat(y)_i)^2, sum^N_(i = 1)(y_i - overline(y))^2) $
В рамках рассматриваемой нами модели нетрудно увидеть, что $overline(y) = display(frac(n_1, N))$, где $n_1$ — число наблюдений, для которых $y_i = 1$.
Если раскрыть в знаменателе скобки, то мы получим
$
  sum^N_(i = 1)(y_i - overline(y))^2 & = sum^N_(i = 1)y^2_i - 2 sum^N_(i = 1)y_i overline(y) + sum^N_(i = 1)overline(y)^2 = \
  & = n_1 - 2 overline(y) sum^N_(i = 1)y_i + N(frac(n_1, N))^2 = \
  & = n_1 - 2 frac(n^2_1, N) + frac(n^2_1, N) = frac(n_0 n_1, N)
$

Подставив это значение в формулу для $R^2$, получим коэффициент детерминации Эфрона #cite(<efron1978regression>)
$ R^2 = 1 - frac(N, n_0 n_1)sum^N_(i = 1)(y_i - hat(y)_i)^2 $

С другой стороны, вспомним, что рассматриваемые модели оцениваются при помощи метода максимального правдоподобия, поэтому хотелось бы
как-то использовать этот факт.
В стандартной линейной регрессии @baguley2012online@pierse2021limited
#eqnum[$ R^2 = 1 - (frac(L_0, L))^(2 slash N) $] <eq-ch5-R2>
где $L$ — максимум функции правдоподобия для модели, а $L_0$ — максимум функции правдоподобия для модели, в которой единственным
регрессором является константа.

Однако в моделях бинарного выбора
$L_0 lt.eq.slant L lt.eq.slant 1$, поэтому $R^2 lt.eq.slant 1 - L^(2 slash N)_0$,
и нам необходимо как-то перемасштабировать эту меру.

Крэгг и Улер #cite(<cragg1970demand>) предложили следующую меру
$ R^2_"pseudo" = frac(L^(2 slash N) - L^(2 slash N)_0, (1 - L^(2 slash N)_0) L^(2 slash N)), $
которая, на самом деле, представляет собой формулу @eq-ch5-R2, деленную на $1 - L^(2 slash N)_0$.

МакФадден #cite(<mcfadden1974measurement>) предложил другую меру
$ R^2_"McFadden" = 1 - frac(ln L, ln L_0) $

Разумеется, что можно использовать информационные критерии
$
   italic("AIC") & = - frac(2 L, N) + frac(2 p, N) \
   italic("BIC") & = - frac(2 L, N) + frac(p ln p, N) \
  italic("HQIC") & = - frac(2 L, N) + frac(2 p ln ln N, N)
$

= Дискретные модели

== Порядковые модели

Предположим теперь, что
$
  y_i = cases(
    1\, & quad y^(*)_i lt.eq.slant gamma_(i\,1)\,,
    dots.v & quad dots.v,
    k\, & quad gamma_(i\,k - 1)< y^(*)_i lt.eq.slant gamma_(i\,k)\,,
    dots.v & quad dots.v,
    K\, & quad gamma_(i\,K - 1)< y^(*)_i
  )
$
где $gamma_(i,1) < ... < gamma_(i,k) < ... < gamma_(i,K - 1)$ — пороговые значения, в общем случае ненаблюдаемые.

Для простоты предположим, что пороговые значения одинаковы для всех наблюдений, то есть $gamma_(i,k) = gamma_k, forall i$.
Для возможности идентификации всех параметров модели обычно также вводится нормализация $Var[u_i] = 1$, а также $gamma_1 = 0$.

В этом случае имеем вероятности
$
  Pr (y_i = 1) & = Pr (y^(*)_i lt.eq.slant 0) = F(- vv(x)^prime_i vv(beta)), \
  Pr (y_i = 2) & = Pr (0 < y^(*)_i lt.eq.slant gamma_2) = F(gamma_2 - vv(x)^prime_i vv(beta)) - F(- vv(x)^prime_i vv(beta)), \
  Pr (y_i = 3) & = Pr (gamma_2< y^(*)_i lt.eq.slant gamma_3) = F(gamma_3 - vv(x)^prime_i vv(beta)) - F(gamma_2 - vv(x)^prime_i vv(beta)), \
  ...
$
на основе которых можно записать функцию правдоподобия
$
  L = product_(y_i = 1) Pr (y_i = 1)... product_(y_i = K) Pr (y_i = K)
$

== Мультиномиальные модели

В предыдущем случае значения $y_i$ были упорядочены, примером могут служить шкалы "от одного до пяти".
Однако бывают ситуации, когда варианты невозможно упорядочить, например, профессии.

Для внесения порядка в эту ситуацию считают, что каждой альтернативе $j$ индивидом $i$ назначается функция полезности $U_(i j)$, равная
$ U_(i j) = vv(x)^prime_(i j)vv(beta) + epsilon_(i j), $
где значения регрессоров могут зависеть не только от наблюдения, но и от альтернативы.
В данной постановке задачи индивид $i$ выбирает альтернативу $j$ с вероятностью
$ Pr (y_i = j) = Pr (U_(i j) = max { U_(i 1),...,U_(i K) }) $

Как отмечает Носко #cite(<nosko2021ekonometrikaa>), в явном виде эту вероятность выразить довольно сложно.
Однако, если ошибки $epsilon_(i j)$ имеют одинаковое логарифмическое распределение Вейбулла (или распределение Гумпеля)
$ F (x) = exp { - e^(- x) }, $
то искомая вероятность будет иметь весьма простой вид
$ Pr (y_i = j) = frac(exp { vv(x)^prime_(i j)vv(beta) }, sum_j exp { vv(x)^prime_(i j)vv(beta) }), $
где для нормализации модели считается, что $vv(x)_(i 1) = vv(0)$.

Если $K = 2$, то модель превращается в логит модель.
#v(2em)

= Урезанные и цензурированные модели

== Цензурированные модели

Рассмотрим задачу, которую решал Тобин #cite(<tobin1958estimation>): он исследовал расходы домохозяйств на покупку автомобиля.
Для некоторых семей эти расходы были равны 0, то есть был отказ от покупки автомобиля.
$ p_i = cases(p^(*)_i &" если семья имеет автомобиль", 0 &" иначе") $
#align(center)[#image("fig/tobit_censored.pdf")]

В этом случае оценки регрессий расходов на доходы будут смещенными и несостоятельными: подобранное значение в модели регрессии показывает цену, которую семья готова заплатить на своем уровне доходов, однако нулевое значение может означать, что на рынке просто _нет_ автомобилей по такой цене.

Предположим, что i-я семья покупает автомобиль по цене $p^(*)_i$ если последняя превышает минимально возможную стоимость $gamma$ автомобиля на рынке (первичном и вторичном), т.е. если $p^(*)_i > gamma$. Для остальных семей известно, что $p^(*)_i <= gamma$. Такие данные называют "цензурированными" на уровне $gamma$, причем наблюдаемой является переменная
$ p^c_i = cases(p^(*)_i &" если " p^(*)_i > gamma, gamma &" если " p^(*)_i <= gamma) $
Данную переменную можно скорректировать на значение $gamma$ и получить
$
  y_t = cases(p^(*)_i - gamma &" если " p^(*)_i > gamma, 0 &" если " p^(*)_i <= gamma) = cases(y^(*)_i &" если " y^(*)_i > 0, 0 &" если " y^(*)_i <= 0)
$
где $y^(*)_i$ --- латентная переменная, для которой можно построить *стандартную тобит-модель*

$
                        vv(y)^(*) & = vv(X) vv(beta) + vv(u) #<eq-ch5-ystar.reg> \
  Ex (y^(*)_i mid(bar.v) vv(x)_i) & = vv(x)^prime_i vv(beta) \
                            y^*_i & tilde.op N (vv(x)^prime_i vv(beta), sigma^2)
$

=== Усеченная регрессия

Предположим, что для расчетов используются только те наблюдения, для которых $y_i > 0$
$ vv(y) = vv(X) vv(beta) + vv(u) $

В такой модели для $w>0$ имеем
$
  Pr (y_i <= w) & = Pr (y^*_i <= w mid(bar.v) y^*_i > 0) = (Pr (0 < y^*_i <= w)) / (Pr (y^*_i > 0)) \
  Pr (0 < y^*_i <= w) & = Pr ((-vv(x)^prime_i vv(beta))/sigma < (y^*_i - vv(x)^prime_i vv(beta))/sigma <= (w - vv(x)^prime_i vv(beta))/sigma) = \
  & = Phi ((w-vv(x)^prime_i vv(beta))/sigma) - Phi ((-vv(x)^prime_i vv(beta))/sigma) \
  Pr (y^*_i > 0) & = 1 - Pr ((y^*_i - vv(x)^prime_i vv(beta))/sigma <= (-vv(x)^prime_i vv(beta))/sigma) = 1 - Phi ((-vv(x)^prime_i vv(beta))/sigma) = Phi ((vv(x)^prime_i vv(beta))/sigma)
$

Отсюда можно найти, что @nosko2021ekonometrikaa @magnus2021ekonometrika
$
  p_y (w) = (dif Pr (y_i <= w)) / (dif w) = 1/sigma frac(phi((w - vv(x)^prime_i vv(beta))/sigma), Phi ((vv(x)^prime_i vv(beta))/sigma))
$
$
  Ex (y_i mid(bar.v) vv(x)_i) & = integral_0^oo w p_y (w) dif w = \
  & = 1/(sigma Phi ((vv(x)^prime_i vv(beta))/sigma)) integral_0^oo w phi((w - vv(x)^prime_i vv(beta))/sigma) dif w
$
#text(size: 10pt)[
  $
    integral_0^oo w phi((w - vv(x)^prime_i vv(beta))/sigma) dif w & = vv(x)^prime_i vv(beta) integral_0^oo phi((w - vv(x)^prime_i vv(beta))/sigma) dif w + sigma integral_0^oo (w - vv(x)^prime_i vv(beta))/sigma phi((w - vv(x)^prime_i vv(beta))/sigma) dif w = \
    & = sigma vv(x)^prime_i vv(beta) integral_(-vv(x)^prime_i vv(beta) slash sigma)^oo phi(t) dif t + sigma^2 integral_(-vv(x)^prime_i vv(beta) slash sigma)^oo t phi(t) dif t = \
    & = sigma vv(x)^prime_i vv(beta) [1 - Phi(-(vv(x)^prime_i vv(beta))/sigma)] + sigma^2/sqrt(2 pi) integral_(-vv(x)^prime_i vv(beta) slash sigma)^oo t exp{-t^2/2} dif t \
    & = vv(x)^prime_i vv(beta) sigma Phi((vv(x)^prime_i vv(beta))/sigma) - sigma^2/(sqrt(2 pi)) lr(exp{-t^2/2}|)^oo_(-vv(x)^prime_i vv(beta) slash sigma) = \
    & = vv(x)^prime_i vv(beta) sigma Phi((vv(x)^prime_i vv(beta))/sigma) + sigma^2/sqrt(2 pi) exp{-(-vv(x)^prime_i vv(beta) slash sigma)^2/2} = \
    & = vv(x)^prime_i vv(beta) sigma Phi((vv(x)^prime_i vv(beta))/sigma) + sigma^2 phi(-(vv(x)^prime_i vv(beta))/sigma)
  $
]
Откуда
$
  Ex (y_i mid(bar.v) vv(x)_i) = vv(x)^prime_i vv(beta) + sigma frac(phi((vv(x)^prime_i vv(beta))/sigma), Phi((vv(x)^prime_i vv(beta))/sigma)) > vv(x)^prime_i vv(beta) #<eq-ch5-EYt>
$

=== Учет нулевых значений

Ту же модель можно оценить и иным способом, учитывая неполные наблюдения.
Объясняемой переменной в рамках данного подхода будет
$
  y_t = cases(
    vv(x)^prime_i vv(beta) & " если " vv(x)^prime_i vv(beta) + u_i > 0,
    0 & " если " vv(x)^prime_i vv(beta) + u_i <= 0,
  )
$

В этом случае имеем, что
$
  Pr (y_i = 0 mid(bar.v) vv(x)_i) & = Pr (u_i <= -vv(x)^prime_i vv(beta)) = Phi (-(vv(x)^prime_i vv(beta))/sigma) = 1 - Phi ((vv(x)^prime_i vv(beta))/sigma) \
  Pr (y_i > 0 mid(bar.v) vv(x)_i) & = Phi ((vv(x)^prime_i vv(beta))/sigma) \
  Pr (y_i <= w mid(bar.v) vv(x)_i) & = Pr ((y_i - vv(x)^prime_i vv(beta))/sigma <= (w - vv(x)^prime_i vv(beta))/sigma) = Phi((w - vv(x)^prime_i vv(beta))/sigma)
$

Отсюда мы можем найти математическое ожидание $y_i$:
$
  Ex (y_i mid(bar.v) vv(x)_i) & = 0 times Pr (y_i = 0 mid(bar.v) vv(x)_i) + integral_0^oo w phi((w - vv(x)^prime_i vv(beta))/sigma) dif w times Pr (y_i > 0 mid(bar.v) vv(x)_i) = \
  & = integral_0^oo w phi((w - vv(x)^prime_i vv(beta))/sigma) dif w times Phi((vv(x)^prime_i vv(beta))/sigma) = \
  & = [vv(x)^prime_i vv(beta) + sigma frac(phi((vv(x)^prime_i vv(beta))/sigma), Phi((vv(x)^prime_i vv(beta))/sigma))] times Phi((vv(x)^prime_i vv(beta))/sigma) #<eq-ch5-EY>
$

Из полученного выражения понятно, что оценка модели
#eqnum[$ vv(y) = vv(x)^prime_i vv(beta) + vv(u) $] <eq-ch5-y.reg>
будет смещена и несостоятельна, так как $Ex (vv(u)) != 0$.

На следующем рисунке это наглядно показано.
Черные и серые закрашенные кружки обозначают наблюдения $(x_i,y^(*)_i)$, незакрашенные кружки показывают цензурированные наблюдения $(x_i,0)$ для $y^(*)_i< 0$.
Черная штриховая линия показывает "истинную" линию регрессии @eq-ch5-ystar.reg, красная линия — подобранную линию регрессии @eq-ch5-y.reg, а синяя — график $Ex (y_i mid(bar.v) vv(x)_i)$, построенный по формуле @eq-ch5-EY.
Видно, что красная линия смещена относительно истинной линии.
Синяя линия, в свою очередь, по мере удаления от области с цензурированными значениями стремится к истинной линии.
#align(center)[#image("fig/tobit.pdf")]

Для получения оценок следует воспользоваться методом максимального правдоподобия.
Функция правдоподобия имеет вид
$
  L= product_(y_i = 0)[ 1 - Phi (frac(vv(x)^prime_i vv(beta), sigma)) ] product_(y_i> 0)frac(1, sqrt(2 pi sigma^2))exp { - frac((y_i - vv(x)^prime_i vv(beta))^2, 2 sigma^2) }
$

Предельный эффект объясняющих факторов будет равен
$
  frac(partial Ex (y_i mid(bar.v) vv(x)_i), partial vv(x)_i) = Phi (frac(vv(x)^prime_i vv(beta), sigma)) vv(beta)
$

= Усеченные модели

Рассмотрим еще один вариант: ситуацию, когда при определенных условиях наблюдения вообще недоступны.
Представим себе обычную модель
$ vv(y) = vv(X) vv(beta) + vv(u), vv(u) tilde.op N (vv(0),sigma^2 vv(I)) $
и пусть нам доступны только те наблюдения, для которых $y_i> a$.
Тогда из @eq-ch5-EYt следует
$
  Ex (y_i mid(bar.v) y_i > a) & = vv(x)^prime_i vv(beta) + sigma frac(
                                  phi (frac(a - vv(x)^prime_i vv(beta), sigma)),
                                  Phi (frac(a - vv(x)^prime_i vv(beta), sigma))
                                ) = \
                              & = vv(x)^prime_i vv(beta) + sigma lambda (frac(a - vv(x)^prime_i vv(beta), sigma))
$
то есть если для данной урезанной выборки оценить модель при помощи МНК, то оценки будут смещены, поскольку мы получим, что $Ex (u_i mid(bar.v) vv(x)_i) = sigma lambda (...)$.
Поэтому состоятельные оценки можно получить при помощи метода максимального правдоподобия
$
  Pr (y_i mid(bar.v) y_i> a) & = 1/sigma frac(phi (frac(y_i - vv(x)^prime_i vv(beta), sigma)), 1 - Phi (frac(a - vv(x)^prime_i vv(beta), sigma))) \
  L & = frac(1, sigma^N)[ 1 - Phi (frac(a - vv(x)^prime_i vv(beta), sigma)) ]^N product_i frac(1, sqrt(2 pi sigma^2))exp { - frac((y_i - vv(x)^prime_i vv(beta))^2, 2 sigma^2) }
$

#bibliography("bibliography.bib")

#show: appendix

= Математическое ожидание усеченного нормального распределения

Пусть $y^* tilde.op N(mu, sigma^2)$.
Тогда случайная величина
$
  y = cases(
    y^* & " если " y_t in [a, b],
    0 & " иначе"
  )
$
будет иметь усеченное нормальное распределение.

Из определения плотности следует, что
$
  f(x) = cases(
    display(1/sigma frac(phi((x-mu)/sigma), Phi((b-mu)/sigma) - Phi((a-mu)/sigma))) & " если " x in [a, b],
    0 & " иначе"
  )
$

Отсюда
$
  Ex y & = integral_(-oo)^oo x f(x) dif x = integral_a^b x f(x) dif x = \
  & = 1/sigma frac(1, Phi((b-mu)/sigma) - Phi((a-mu)/sigma)) underbrace(integral_a^b x phi((x-mu)/sigma) dif x, Q) = bigast
$
Рассмотрим интеграл $Q$
$
  Q & = integral_a^b (x #text(fill: red)[$plus.minus mu$]) phi((x-mu)/sigma) dif x = \
  & = mu integral_a^b phi((x-mu)/sigma) dif x + #text(fill: blue)[$sigma$] integral_a^b (x-mu)/#text(fill: blue)[$sigma$] phi((x-mu)/sigma) dif x = \
  & = mat(delim: bar.v, t = display((x-mu)/sigma), x = mu + sigma t, dif x = sigma dif t) = \
  & = sigma mu integral_((a-mu)/sigma)^((b-mu)/sigma) phi(t) dif t + sigma^2 integral_((a-mu)/sigma)^((b-mu)/sigma) t phi(t) dif t = \
  & = sigma mu integral_((a-mu)/sigma)^((b-mu)/sigma) dif Phi(t) - sigma^2/sqrt(2 pi) integral_((a-mu)/sigma)^((b-mu)/sigma) dif exp{-t^2/2} = \
  & = sigma mu [Phi((b-mu)/sigma) - Phi((a-mu)/sigma)] + sigma^2 [exp{-((b-mu)\/sigma)^2/2}/sqrt(2 pi) - exp{-((a-mu)\/sigma)^2/2}/sqrt(2 pi)] = \
  & = sigma mu [Phi((b-mu)/sigma) - Phi((a-mu)/sigma)] + sigma^2 [phi((b-mu)/sigma) - phi((a-mu)/sigma)]
$
Здесь мы воспользовались тем, что $dif e^(-t^2/2) = -t e^(-t^2/2) dif t$.
Откуда получим
$
  bigast = mu + sigma frac(phi((b-mu)/sigma) - phi((a-mu)/sigma), Phi((b-mu)/sigma) - Phi((a-mu)/sigma))
$

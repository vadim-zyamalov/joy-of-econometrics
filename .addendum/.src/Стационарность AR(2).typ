#import "@local/trinkets:0.1.0": *
#import "@local/neat-document:0.1.0": *

#set document(title: [Стационарность AR(2)])
#show: neat-document

Характеристическое уравнение для данного процесса имеет вид
$ lambda^2 - a_1 lambda - a_2 = 0 $

Вид решения для однородной части зависит от дискриминанта
$ D = a^2_1 + 4 a_2 $

/ $D > 0$: корни $lambda_(1,2) in RR$, причем $lambda_(1)!= lambda_2$.
  Решение однородного уравнения имеет вид
  $ y_t = C_1 lambda^t_1 + C_2 lambda^t_2 $

  Очевидно, что оно стабильно при $abs(lambda_(1,2)) < 1$.

/ $D = 0$: корни $lambda_(1,2) in RR$, причем $lambda_1 = lambda_2$.
  Решение однородного уравнения имеет вид
  $ y_t = C_1 lambda^t_1 + C_2 t lambda^t_1 $

  Оно стабильно при $abs(lambda_1) < 1$.

/ $D < 0$: корни $lambda_(1,2) in CC$, причем $lambda_1 = lambda^*_2$.

  Комплексные числа имеют три эквивалентные записи
  $ z = a + i b = r e^(i theta)= r[ cos (theta) + i sin (theta) ] $
  где $r = sqrt(a^2 + b^2)$ — модуль числа, $theta = arctan b/a$ -- аргумент числа.
  Тогда комплексно сопряженные корни имеют вид $lambda_(1,2)= r e^(pm i theta)$.

  Отсюда получим, что
  $
    y_t & = C_1 lambda^t_1 + C_2 lambda^t_2 = \
        & = C_1 r^t e^(i theta t)+ C_2 r^t e^(- i theta t)= \
        & = C_1 r^(t)[ cos (theta t) + i sin (theta t) ] + C_2 r^(t)[ cos (theta t) - i sin (theta t) ] = \
        & = (C_1 + C_2) r^t cos (theta t) + (C_1 - C_2) i sin (theta t)
  $

  Вспомним, что $Y_t$ должно быть _вещественным_!
  Для этого нужно, чтобы $C_1 + C_2 in RR$ и $C_1 - C_2 in CC$.
  Это возможно только если $C_1 = C^*_2$. То есть $C_(1,2) = v e^(pm i xi)$.

  Получим
  $
    y_t & = v r^t e^(i(theta t + xi))+ v r^t e^(- i(theta t + xi))= \
        & = v r^t (e^(i(theta t + xi))+ e^(- i(theta t + xi))) = \
        & = 2 v r^t cos (theta t + xi)
  $

  Последнее следует из
  $
    e^(i x)+ e^(- i x) & = cos (x) + i sin (x) + cos (- x) + i sin (- x) = \
                       & = cos (x) + i sin (x) + cos (x) - i sin (x) = \
                       & = 2 cos (x)
  $

  Значение $v$ зависит от начальных условий и неоднородного члена уравнения.
  Очевидно, что решение стабильно, если $r < 1$.

#align(center)[#image("fig/AR2.sols.pdf")]

== Коэффициенты процесса

Если $D = a^2_1 + 4 a_2 > 0$, то есть $a_2 > - 1/4 a^2_1$, то
$
  abs(frac(a_1 pm sqrt(a^2_1+ 4 a_2), 2)) < 1 \
  -1 < frac(a_1 - sqrt(a^2_1+ 4 a_2), 2) < frac(a_1 + sqrt(a^2_1+ 4 a_2), 2) < 1
$

Из первого неравенства получим
$
                  - 2 & < a_1 - sqrt(a^2_1 + 4 a_2) \
  sqrt(a^2_1 + 4 a_2) & < a_1 + 2
$
откуда следует, что
$ cases(a_1 + 2 > 0, a^2_1 + 4 a_2 < a^2_1 + 4 a_1 + 4) => cases(a_1 > -2, a_2 < a_1 + 1) $

Из третьего неравенства получим
$
                    2 & > a_1 + sqrt(a^2_1 + 4 a_2) \
  sqrt(a^2_1 + 4 a_2) & < 2 - a_1
$
откуда следует, что
$ cases(2 - a_1 > 0, a^2_1 + 4 a_2 < a^2_1 - 4 a_1 + 4) => cases(a_1 < 2, a_2 < - a_1 + 1) $

#align(center)[#image("fig/AR2.sols.g.pdf")]

Если $D = 0$, то есть $a_2 = - 1/4 a^2_1$, то
$
  abs(frac(a_1, 2)) < 1 \
  - 2 < a_1 < 2
$

#align(center)[#image("fig/AR2.sols.e.pdf")]

Если же $D < 0$, то есть $a_2 < -1/4 a^2_1$, то
$
  lambda_(1,2) & = frac(a_1, 2) pm i frac(sqrt(-4 a_2 - a^2_1), 2) \
             r & = sqrt(Re^2 lambda + Im^2 lambda) \
             r & = sqrt(frac(a^2_1, 4)- frac(4 a_2 + a^2_1, 4))= sqrt(-a_2) \
             0 & < sqrt(-a_2) < 1 \
            -1 & < a_2 < 0
$
#align(center)[#image("fig/AR2.sols.l.pdf")]


Все это можно подытожить @fig-p4.6.
#figure(
  [#image("fig/AR2.sols.all.pdf", width: 100%)],
  caption: [Область стационарности процесса AR(2)],
) <fig-p4.6>

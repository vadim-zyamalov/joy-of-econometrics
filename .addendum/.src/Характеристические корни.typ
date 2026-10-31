#import "@local/trinkets:0.1.0": *
#import "@local/neat-document:0.1.0": *

#set document(title: [Корни характеристического уравнения и разложение лагового полинома AR(p)])
#show: neat-document

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

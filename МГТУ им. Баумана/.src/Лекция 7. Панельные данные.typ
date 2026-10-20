#import "@preview/theorion:0.6.0": *
#import "@preview/tschich:0.2.0": *
#import "@local/neat-document:0.1.0": *
#import "@local/trinkets:0.1.0": *
#set document(title: [Панельные данные])
#let margins = tschich-var(210mm, 297mm, 1 / 10)
#let font-size = 14pt
#show: neat-document.with(margin: margins, leading: 1em, fonts: (size: font-size))
#show: show-theorion
#show: trinkets

Ранее мы с вами рассматривали
- межобъектные данные, представляющие собой набор наблюдений для разных объектов в один момент времени#footnote[Либо для разных моментов, но при условии несущественности различий во времени.], и
- временные ряды, представляющие собой наблюдения для одного объекта в последовательные моменты времени.

Панельные данные являются объединением этих двух концепций: набор данных представляет собой наблюдения для большого числа объектов в последовательные моменты времени.

= Основные модели

== Модель сквозной регрессии (пула)

Простейшей моделью панельных данных является модель пула, которая в покомпонентной записи имеет вид
$ y_(i t) = Vect(x)^(prime)_(i t)Vect(beta) + mu + epsilon_(i j), $
а в матричном
$
  underbrace(
    mat(
      delim: "(", Vect(y)_(1);
      Vect(y)_(2);
      dots.v;
      Vect(y)_(N)
    ), Vect(y)
  ) = underbrace(
    mat(
      delim: "(", Vect(X)_(1);
      Vect(X)_(2);
      dots.v;
      Vect(X)_(N)
    ), Vect(X)
  ) Vect(beta) + mu Vect(1)_(N T) + underbrace(
    mat(
      delim: "(", Vect(epsilon)_(1);
      Vect(epsilon)_(2);
      dots.v;
      Vect(epsilon)_(N)
    ), Vect(epsilon)
  ) = underbrace(
    mat(
      delim: "(", Vect(1)_T, Vect(X)_(1);
      Vect(1)_T, Vect(X)_(2);
      dots.v, dots.v;
      Vect(1)_T, Vect(X)_(N)
    ), tilde(Vect(X))
  ) mat(
    delim: "(", mu ;
    Vect(beta)
  ) + Vect(epsilon),
$
где $Vect(y)_i = mat(delim: "(", y_(i 1), y_(i 2), ..., y_(i T))^(prime)$,
$Vect(X)_i = mat(delim: "(", Vect(x)_(i 1), Vect(x)_(i 2), ..., Vect(x)_(i T))^(prime)$.

Оценка коэффициентов сквозной регрессии получается при помощи МНК по обычной формуле
$
  mat(
    delim: "(", hat(mu);
    hat(Vect(beta))
  ) = (tilde(Vect(X))^(prime)tilde(Vect(X)))^(- 1)tilde(Vect(X))^(prime)Vect(y)
$

Данная модель, будучи самой простой, является самой ограниченной, так как она предполагает, что все объекты в выборке являются гомогенными, не обладающими индивидуальными особенностями.
#align(center)[#image("fig/panel-pool.pdf", width: 10cm)]

В случае, когда объекты являются гетерогенными, то есть обладающими набором индивидуальных особенностей, оценка модели пула, в общем случае, будет смещенной.
#align(center)[#image("fig/panel-FE.pdf", width: 10cm)]

== Модель с фиксированными эффектами

Попытаемся учесть данную гетерогенность, разрешив свободному члену принимать различные значения для разных объектов, что в покомпонентной записи будет иметь вид
$
  y_(i t) = Vect(x)^(prime)_(i t)Vect(beta) + mu_i + u_(i j),
  #<eq-ch7-FE>
$
а в матричной
$
  Vect(y) & = Vect(X) Vect(beta) + overbrace(
              mat(
                delim: "(", Vect(1)_T, Vect(0)_T, dots.h.c, Vect(0)_(T);
                Vect(0)_T, Vect(1)_T, dots.h.c, Vect(0)_(T);
                dots.v, dots.v, dots.down, dots.v;
                Vect(0)_T, Vect(0)_T, dots.h.c, Vect(1)_(T)
              ), N
            ) Vect(mu) + Vect(epsilon) = \
          & = Vect(X) Vect(beta) + (Vect(I)_N times.o Vect(1)_T) Vect(mu) + Vect(epsilon)
$
где $Vect(1)_T = underbrace(mat(delim: "(", 1, ..., 1)^(prime), T)$,
$Vect(0)_T = underbrace(mat(delim: "(", 0, ..., 0)^(prime), T)$,
$Vect(mu) = mat(delim: "(", mu_1, ..., mu_N)^(prime)$.

Условимся при этом, что $Vect(X)$ не включает в себя свободный член, дабы избежать ловушки дамми переменных#footnote[Разумеется, можно включить константу в $Vect(X)$ и наложить ограничения на $Vect(mu)$, однако, как мы увидим далее, это приведет к проблемам с вычислениями.].
Данная модель оценивается при помощи обычного МНК и в литературе называется LSDV—МНК с фиктивными переменными.

Очевидным слабым местом данной модели является потенциально большое количество оцениваемых параметров, что приводит к необходимости обращения матрицы очень большой размерности.
С другой стороны можно воспользоваться подходом частичной регрессии для получения оценок коэффициентов $Vect(beta)$ @ratnikova2006vvedenieb@pierse2021panel:
$
  hat(Vect(beta)) = (Vect(X)^(prime)Vect(M)_(mu)Vect(X))^(- 1)Vect(X)^(prime)Vect(M)_(mu)Vect(y),
  #<eq-ch7-FE.estimator>
$
где
$
  Vect(M)_(mu) & = Vect(I)_(N T) - (Vect(I)_N times.o Vect(1)_T) [ (Vect(I)_N times.o Vect(1)_T)^(prime)(Vect(I)_N times.o Vect(1)_T) ]^(- 1)(Vect(I)_N times.o Vect(1)_T)^(prime) = \
  & = Vect(I)_(N T) - (Vect(I)_N times.o Vect(1)_T) [ (Vect(I)_N times.o Vect(1)^(prime)_T) (Vect(I)_N times.o Vect(1)_T) ]^(- 1)(Vect(I)_N times.o Vect(1)_T)^(prime) = \
  & = Vect(I)_(N T) - (Vect(I)_N times.o Vect(1)_T) [ Vect(I)_N times.o Vect(1)^(prime)_T Vect(1)_T ]^(- 1)(Vect(I)_N times.o Vect(1)_T)^(prime) = \
  & = Vect(I)_(N T) - 1/T (Vect(I)_N times.o Vect(1)_T) (Vect(I)_N times.o Vect(1)_T)^(prime) = \
  & = Vect(I)_(N T) - 1/T (Vect(I)_N times.o Vect(1)_T Vect(1)^(prime)_T) = \
  & = Vect(I)_N times.o(Vect(I)_T - frac(Vect(J)_T, T)) = Vect(W)
$
где $Vect(J)_T$ — матрица $T times T$, полностью заполненная единицами.
При выводе этой формулы использовались свойства произведения Кронеккера.

Разберем оператор $Vect(W)$ подробнее:
+ $frac(Vect(J)_T, T)$, будучи примененной к подходящему по размеру вектору $Vect(v)$, дает нам вектор аналогичной длины, заполненный выборочными средними значений исходного вектора, то есть $macron(v) Vect(1)_T$.
+ $Vect(I)_T - frac(Vect(J)_T, T)$, будучи примененной к подходящему по размеру вектору $Vect(v)$, дает нам вектор аналогичной длины, заполненный отклонениями значений элементов вектора от их выборочного среднего.
+ Кронеккерово произведение создает блочную матрицу, применение которой, например, к $Vect(y)$ дает нам вектор центрированных значений, причем для каждого объекта значения центрируются относительного своего выборочного среднего.

Нетрудно увидеть, что данная процедура аналогична таковой, приведенной, например, в @nosko2021ekonometrikaa.
Полученная оценка $hat(Vect(beta))$ будет несмещенной и состоятельной во времени; ее называют "внутри"-оценкой, или _within_-оценкой#footnote[Поэтому мы и назвали оператор $Vect(W)$.], или FE-оценкой.

Индивидуальные эффекты оцениваются по формуле
$ mu_i = overline(Vect(y))_i - overline(Vect(X))_i hat(Vect(beta))_(F E) $

Для дальнейшего нам будет удобно ввести оператор
$ Vect(B) = Vect(I)_N times.o frac(Vect(J)_T, T), $
который, будучи подставлен в @eq-ch7-FE.estimator даст нам несмещенную и состоятельную во времени "между"-оценку (или _between_-оценку).
Нетрудно убедиться в том, что операторы $Vect(W)$ и $Vect(B)$ симметричны и идемпотентны; более того, они являются ортогональными
$ Vect(W) Vect(B) = Vect(B) Vect(W) = Vect(0), $
а также $Vect(W) + Vect(B) = Vect(I)_(N T)$.

== Модель со случайными эффектами

Предположим теперь, что индивидуальные эффекты в @eq-ch7-FE теперь не константы, а случайные величины, некоррелированные с регрессорами и ошибками $Vect(epsilon)$.

Перепишем модель в виде
$
  Vect(y) & = Vect(X) Vect(beta) + Vect(u) \
  u_(i t) & = mu_i + epsilon_(i t)
$

Исходя из предположений относительно ошибок и случайных эффектов
$
  Ex (Vect(epsilon) Vect(epsilon)^(prime)) & = sigma^2_(epsilon)Vect(I)_(N T), \
            Ex (Vect(mu) Vect(mu)^(prime)) & = sigma^2_(mu)Vect(I)_N, \
       Ex (Vect(mu) Vect(epsilon)^(prime)) & = Vect(0),
$
ковариационная матрица для $Vect(u)_i$ имеет вид
$
  Vect(Sigma) = Ex (Vect(u)_i Vect(u)^(prime)_i) = mat(
    delim: "(", sigma^2_(mu) + sigma^2_(epsilon), sigma^2_(mu), dots.h.c, sigma^2_(mu);
    sigma^2_(mu), sigma^2_(mu) + sigma^2_(epsilon), dots.h.c, sigma^2_(mu);
    dots.v, dots.v, dots.down, dots.v;
    sigma^2_(mu), sigma^2_(mu), dots.h.c, sigma^2_(mu) + sigma^2_(epsilon)
  ) = sigma^2_(mu)Vect(J)_T + sigma^2_(epsilon)Vect(I)_T,
$
а ковариационная матрица $Vect(u)$
$
  Vect(Omega) = mat(
    delim: "(", Vect(Sigma), Vect(0), dots.h.c, Vect(0);
    Vect(0), Vect(Sigma), dots.h.c, Vect(0);
    dots.v, dots.v, dots.down, dots.v;
    Vect(0), Vect(0), dots.h.c, Vect(Sigma)
  ) = Vect(I)_N times.o Vect(Sigma)
$
Данные матрицы можно переписать через введенные выше операторы:
$
  Vect(Sigma) & = sigma^2_(mu)Vect(J)_T + sigma^2_(epsilon)Vect(I)_T = \
  & = T sigma^2_(mu)frac(Vect(J)_T, T) + sigma^2_(epsilon)(Vect(I)_T + frac(Vect(J)_T, T) - frac(Vect(J)_T, T)) = \
  & = (T sigma^2_(mu) + sigma^2_(epsilon)) frac(Vect(J)_T, T) + sigma^2_(epsilon)(Vect(I)_T - frac(Vect(J)_T, T)), \
  Vect(Omega) & = Vect(I)_N times.o Vect(Sigma) = \
  & = (T sigma^2_(mu) + sigma^2_(epsilon)) (Vect(I)_N times.o frac(Vect(J)_T, T)) + sigma^2_(epsilon)[ Vect(I)_n times.o(Vect(I)_T - frac(Vect(J)_T, T)) ] = \
  & = (T sigma^2_(mu) + sigma^2_(epsilon)) Vect(B) + sigma^2_(epsilon)Vect(W) = \
  & = sigma^2_(epsilon)(Vect(W) + frac(1, frac(sigma^2_(epsilon), T sigma^2_(mu) + sigma^2_(epsilon)))Vect(B)) = sigma^2_(epsilon)(Vect(W) + frac(1, theta^2)Vect(B)), \
  Vect(Omega)^(- 1) & = frac(1, sigma^2_(epsilon))(Vect(W) + theta^2 Vect(B))
  #<eq-ch7-omega.inv>
$
Данные формулы показывают, что ошибки модели со случайными эффектами
автокоррелированы, то есть FE-МНК оценки будут неэффективными; для
эффективных оценок следует воспользоваться обобщенным МНК
$
  hat(Vect(beta))_(R E) = (Vect(X)^(prime)Vect(Omega)^(- 1)Vect(X))^(- 1)Vect(X)^(prime)Vect(Omega)^(- 1)Vect(y),
  #<eq-ch7-RE>
$
заменив все неизвестные параметры генеральной совокупности выборочными
аналогами.
Кстати, из формул @eq-ch7-omega.inv и @eq-ch7-RE можно получить следующее представление
$
  hat(Vect(beta))_(R E) & = [ Vect(X)^(prime)(Vect(W) + theta^2 Vect(B)) Vect(X) ]^(- 1)Vect(X)^(prime)(Vect(W) + theta^2 Vect(B)) Vect(y) = \
  & = (Vect(X)^(prime)Vect(W) Vect(X) + theta^2 Vect(X)^(prime)Vect(B) Vect(X))^(- 1)Vect(X)^(prime)Vect(W) Vect(y) + (Vect(X)^(prime)Vect(W) Vect(X) + theta^2 Vect(X)^(prime)Vect(B) Vect(X))^(- 1)theta^2 Vect(X)^(prime)Vect(B) Vect(y) = \
  & = underbrace((Vect(X)^(prime)Vect(W) Vect(X) + theta^2 Vect(X)^(prime)Vect(B) Vect(X))^(- 1)(Vect(X)^(prime)Vect(W) Vect(X)), Vect(omega)) (Vect(X)^(prime)Vect(W) Vect(X))^(- 1)Vect(X)^(prime)Vect(W) Vect(y) + zws \
  & wide + underbrace((Vect(X)^(prime)Vect(W) Vect(X) + theta^2 Vect(X)^(prime)Vect(B) Vect(X))^(- 1)theta^2 (Vect(X)^(prime)Vect(B) Vect(X)), Vect(lambda)) (Vect(X)^(prime)Vect(B) Vect(X))^(- 1)Vect(X)^(prime)Vect(B) Vect(y) = \
  & = Vect(omega) hat(Vect(beta))_W + Vect(lambda) hat(Vect(beta))_B,
  #<eq-ch7-decomp>
$
причем
$
  Vect(omega) + Vect(lambda) & = (Vect(X)^(prime)Vect(W) Vect(X) + theta^2 Vect(X)^(prime)Vect(B) Vect(X))^(- 1)(Vect(X)^(prime)Vect(W) Vect(X)) + (Vect(X)^(prime)Vect(W) Vect(X) + theta^2 Vect(X)^(prime)Vect(B) Vect(X))^(- 1)theta^2 (Vect(X)^(prime)Vect(B) Vect(X)) = \
  & = (Vect(X)^(prime)Vect(W) Vect(X) + theta^2 Vect(X)^(prime)Vect(B) Vect(X))^(- 1)(Vect(X)^(prime)Vect(W) Vect(X) + theta^2 Vect(X)^(prime)Vect(B) Vect(X)) = Vect(I)
$
То есть ОМНК оценка коэффициентов является взвешенным средним внутри- и между-оценок.

=== Интерпретация параметра тета

Дисперсия within-компоненты ошибки $u_(i t)$ равна
$
  op("Var")_W u_(i t) = op("Var") (u_(i t) - macron(u)_i) = op("Var") (epsilon_(i t) - overline(epsilon)_i) = sigma^2_(epsilon)(1 - 1/T)
$
Дисперсия between-компоненты ошибки $u_(i t)$ равна
$
  op("Var")_B u_(i t) = op("Var") macron(u)_i = op("Var") mu_i + op("Var") overline(epsilon)_i = sigma^2_(mu) + frac(sigma^2_(epsilon), T)
  #<eq-ch7-Var.B>
$
Отсюда
$ theta^2 = frac(1, T - 1)frac(op("Var")_W u_(i t), op("Var")_B u_(i t)), $
то есть $theta^2$ с точностью до множителя $display(frac(1, T - 1))$ является отношением дисперсий внутри- и между-компонент дисперсии ошибок.

=== Оценивание параметра тета

Если взять остатки между-регрессии @ratnikova2006vvedenieb
$
  hat(Vect(u))_B & = Vect(B) Vect(y) - Vect(B) Vect(X) hat(Vect(beta))_B = \
  & = Vect(B) Vect(y) - Vect(B) Vect(X) (Vect(X)^(prime)Vect(B) Vect(X))^(- 1)Vect(X)^(prime)Vect(B) Vect(y) = \
  & = [ Vect(I) - Vect(B) Vect(X) (Vect(X)^(prime)Vect(B) Vect(X))^(- 1)Vect(X)^(prime) ] Vect(B) (Vect(X) Vect(beta) + Vect(u)) = \
  & = [ Vect(I) - Vect(B) Vect(X) (Vect(X)^(prime)Vect(B) Vect(X))^(- 1)Vect(X)^(prime) ] Vect(B) Vect(u) = \
  & = Vect(B) Vect(u) - Vect(B) Vect(X) (Vect(X)^(prime)Vect(B) Vect(X))^(- 1)Vect(X)^(prime)Vect(B) Vect(u) = \
  & = [ Vect(I) - underbrace(Vect(B) Vect(X) (Vect(X)^(prime)Vect(B) Vect(X))^(- 1)Vect(X)^(prime)Vect(B), Vect(P)_(B X)) ] Vect(B) Vect(u)
$
$
  Ex [ hat(Vect(u))^(prime)_B hat(Vect(u))_B ] & = Ex tr [ Vect(u)^(prime)Vect(B) (Vect(I) - Vect(P)_(B X)) Vect(B) Vect(u) ] = \
  & = tr [ Vect(B) (Vect(I) - Vect(P)_(B X)) Vect(B) Ex (Vect(u) Vect(u)^(prime)) ] = \
  & = tr [ Vect(B) (Vect(I) - Vect(P)_(B X)) Vect(B) Vect(Omega) ] = \
  & = tr [ Vect(B) (Vect(I) - Vect(P)_(B X)) Vect(B) sigma^2_(epsilon)(Vect(W) + frac(1, theta^2)Vect(B)) ] = \
  & = tr [ frac(sigma^2_(epsilon), theta^2)Vect(B) (Vect(I) - Vect(P)_(B X)) Vect(B) ] = frac(sigma^2_(epsilon), theta^2)tr [ (Vect(I) - Vect(P)_(B X)) Vect(B) ] = \
  & = frac(sigma^2_(epsilon), theta^2)op("rank") (Vect(B) - Vect(P)_(B X)) = (sigma^2_(epsilon) + T sigma^2_(mu)) (N - K)
$
Здесь мы пользуемся свойством следа матрицы $tr (A B) = tr (B A)$ и свойством $B = B B$.

То есть
$
  Ex [ display(frac(hat(Vect(u))_(B)'hat(Vect(u))_B, N - K)) ] & = sigma^2_B = sigma^2_(epsilon) + T sigma^2_(mu) \
  hat(sigma)^2_B & = display(frac(hat(Vect(u))^(prime)_B hat(Vect(u))_B, N - K))
$

Данная формула отличается от @eq-ch7-Var.B на множитель $T$, так как она рассчитывалась для вектора $hat(Vect(u))_B$ длиной $N T$, а @eq-ch7-Var.B — для вектора средних длиной $N$.

Аналогично
$ hat(sigma)^2_W = hat(sigma)^2_(epsilon) = frac(hat(Vect(u))_(W)'hat(Vect(u))_W, N T - N - K) $

В качестве оценки $theta^2$ можно использовать отношение
$ hat(theta)^2 = frac(hat(sigma)^2_W, hat(sigma)^2_B), $
поскольку
$
  plim_(N -> infinity)hat(theta)^2 = frac(plim_(N -> infinity)hat(sigma)^2_W, plim_(N -> infinity)hat(sigma)^2_B) = theta^(2)
$

= Тестирование

== Тест на правомочность использования модели пула

В данном тесте проводится сравнение оценок модели пула
$ Vect(y)_i = Vect(X)_i Vect(beta) + Vect(epsilon)_i $
против более общей модели
$ Vect(y)_i = Vect(X)_i Vect(beta)_i + Vect(epsilon)_i $

Нулевой гипотезой теста является $H_(0): Vect(beta)_i = Vect(beta), forall i$ @pierse2021panel.
При справедливости нулевой гипотезы $Vect(epsilon) tilde.op N (Vect(0),sigma^2_(epsilon)Vect(I)_(N T))$,
и для тестирования гипотезы может использоваться F-статистика
$ frac(N(T - K), (N - 1) K)frac(S S R_r - S S R_u, S S R_u)tilde.op F ((N - 1) K,N(T - K)), $
где $S S R_r$ — сумма квадратов остатков регрессии пула, а $S S R_u$
— сумма квадратов остатков общей модели.
Аналогичное тестирование можно провести с использованием введенной
выше оценки $hat(theta)$ @ratnikova2006vvedenieb: если
ошибка $Vect(u)$ имеет нормальное распределение, то
#math.equation(block: true, numbering: none)[
  $
    hat(sigma)^2_W & tilde frac(sigma^2_(epsilon), N T - N - K)chi^2_(N T - N - K), \
    hat(sigma)^2_B & tilde frac(sigma^2_(epsilon) + T sigma^2_(mu), N - K)chi^2_(N - K), \
    frac(hat(sigma)^2_(B)/ (sigma^2_(epsilon) + T sigma^2_(mu)), hat(sigma)^2_(W)/ sigma^2_(epsilon)) & tilde frac(chi^2_(N - K)/ (N - K), chi^2_(N T - N - K)/ (N T - N - K)) = F(N - K,N T - N - K),
  $
]
то есть при справедливости нулевой гипотезы $sigma^2_(mu) = 0$
$ display(frac(1, hat(theta)^2)) = frac(hat(sigma)^2_W, hat(sigma)^2_B)tilde.op F (N - K,N T - N - K) $

== Тест Хаусмана

Рассмотрим следующие гипотезы
$
  H_(0): & op("Corr") (mu_i,Vect(X)_i) = 0 \
  H_(a): & op("Corr") (mu_i,Vect(X)_i) != 0
$

Если справедлива нулевая гипотеза, то $hat(Vect(beta))_(R E)$ — состоятельная и асимптотически эффективная, а если альтернативная, то она — асимптотически смещенная; $hat(Vect(beta))_(F E)$ является состоятельной при обеих альтернативах.

В работе Мундлака @mundlak1978pooling предлагается рассматривать и FE, и RE модель как модели со случайным эффектом.
Однако в случае FE модели данный эффект коррелирован с регрессорами, что делает ОМНК оценки несостоятельными; при этом внутри-оценка, являющаяся оценкой, _условной_ относительно реализации случайных эффектов, остается состоятельной.

Введем три величины @ratnikova2006vvedenieb
$
  q_1 & = hat(Vect(beta))_W - hat(Vect(beta))_(R E), \
  q_2 & = hat(Vect(beta))_W - hat(Vect(beta))_B, \
  q_3 & = hat(Vect(beta))_(R E) - hat(Vect(beta))_(B)
$

С учетом разложения @eq-ch7-decomp
$
  q_1 & = (Vect(I) - Vect(omega)) (hat(Vect(beta))_W - hat(Vect(beta))_B) = (Vect(I) - Vect(omega)) q_2, \
  q_3 & = Vect(omega) (hat(Vect(beta))_W - hat(Vect(beta))_B) = Vect(omega) q_(2)
$

В силу данной пропорциональности статистика
$ H= q^(prime)(Var q)^(-1)q tilde chi^2_K $
будет иметь одно и то же значение для всех трех значений $q$.

Проще всего использовать внутри- и между-оценки, так как для них, в силу ортогональности $Vect(W)$ и $Vect(B)$ @ratnikova2006vvedenieb @pierse2021panel
$
  op("Var") (hat(Vect(beta))_W - hat(Vect(beta))_B) & = op("Var") hat(Vect(beta))_W + op("Var") hat(Vect(beta))_B \
  op("Var") hat(Vect(beta))_W & = sigma^2_W (Vect(X)^(prime)Vect(W) Vect(X))^(- 1) = sigma^2_(epsilon)(Vect(X)^(prime)Vect(W) Vect(X))^(- 1) \
  op("Var") hat(Vect(beta))_B & = sigma^2_B (Vect(X)^(prime)Vect(B) Vect(X))^(- 1) = (sigma^2_(epsilon) + T sigma^2_(mu)) (Vect(X)^(prime)Vect(B) Vect(X))^(- 1)
$

В оригинальной статье Хаусмана @hausman1978specification приводится следующая формула, основанная на $q_1$:
$
  H= (hat(Vect(beta))_W - hat(Vect(beta))_(R E))^(prime)(op("Var") hat(Vect(beta))_W - op("Var") hat(Vect(beta))_(R E))^(- 1)(hat(Vect(beta))_W - hat(Vect(beta))_(R E))
$

= Динамические модели

Предположим теперь, что перед нами стоит задача оценить следующую модель
$ y_(i t) = Vect(x)^(prime)_(i t)Vect(beta) + gamma y_(i,t - 1) + alpha_i + epsilon_(i t), $
где $epsilon_(i t)tilde #text[iid] (0,sigma^2_(epsilon))$.

Добавление динамики в модель существенно влияет на интерпретацию коэффициентов, поскольку текущее наблюдение определяется не только текущими значениями регрессоров, но и всей предыдущей историей значений @ratnikova2006vvedenie.
Более того, лагированная переменная коррелирует со случайной ошибкой.

Рассмотрим случай детерминированных $alpha_i$; обычная внутри-оценка на конечных выборках будет смещенной и несостоятельной, так как из-за использования усредненных во времени значений
$
  y_(i t) - macron(y)_i = (Vect(x)_(i t) - overline(Vect(x))_i)^(prime)Vect(beta) + gamma (y_(i,t - 1) - macron(y)_i) + (epsilon_(i t) - overline(epsilon)_t)
$
регрессор $y_(i,t - 1) - macron(y)_i$ и ошибка $epsilon_(i t) - overline(epsilon)_t$ коррелированы.
Более того, оценки исходной регрессии также будут смещенными из-за корреляции между $y_(i,t - 1)$ и $alpha_i$.

@anderson1981estimation предложили устранять индивидуальные эффекты при помощи взятия первых разностей
$
  y_(i t) - y_(i,t - 1)eq.colon Delta y_(i t) = Delta Vect(x)^(prime)_(i t)Vect(beta) + gamma Delta y_(t,t - 1) + Delta epsilon_(i,t)
$

Прямая оценка данной модели также будет смещенной; однако в данном случае можно воспользоваться IV, взяв $y_(i,t - 2)$ в качестве инструмента
для $Delta y_(t,t - 1)$: действительно, $y_(i,t - 2)$ коррелирует с $Delta y_(t,t - 1)$, но не коррелирует с $Delta epsilon_(i,t)$.

== Модель _без_ экзогенных регрессоров

Андерсон и Сяо предложили два варианта состоятельных оценок для $gamma$ в модели без экзогенных регрессоров $Vect(X)$
$ y_(i t) = gamma y_(i,t - 1) + alpha_i + epsilon_(i t) $
с необходимыми условиями состоятельности
$
  plim_(N -> infinity,T -> infinity)frac(1, N (T - 1))sum^N_(i = 1)sum^T_(t = 2) Delta epsilon_(i,t)y_(i,t - 2) = 0
  #<eq-ch7-AB.cond1>
$
$
  plim_(N -> infinity,T -> infinity)frac(1, N (T - 2))sum^N_(i = 1)sum^T_(t = 3) Delta epsilon_(i,t)Delta y_(i,t - 2) = 0
  #<eq-ch7-AB.cond2>
$
Можно заметить, что данные равенства можно переписать как
$
  EE (Delta epsilon_(i,t)y_(i,t - 2)) = 0, \
  EE (Delta epsilon_(i,t)Delta y_(i,t - 2)) = 0,
$
то есть их можно рассматривать как моментные условия в рамках _метода моментов_.

@arellano1991tests предложили расширить количество инструментов другими лагами $y_(i t)$.
Например, для $T = 4$ @ratnikova2006vvedenie
$
  t & = 2 : & Ex (Delta epsilon_(i,2)y_(i 0)) & = 0, \
  t & = 3 : & Ex (Delta epsilon_(i,3)y_(i 0)) & = 0, \
    &       & Ex (Delta epsilon_(i,3)y_(i 1)) & = 0, \
  t & = 4 : & Ex (Delta epsilon_(i,4)y_(i 0)) & = 0, \
    &       & Ex (Delta epsilon_(i,4)y_(i 1)) & = 0, \
    &       & Ex (Delta epsilon_(i,4)y_(i 2)) & = 0
$

Если ввести обозначения
$
  Delta Vect(epsilon)_i & = mat(
                            delim: "(", Delta epsilon_(i 2);
                            Delta epsilon_(i 2);
                            dots.v;
                            Delta epsilon_(i T)
                          ) = Delta Vect(y)_i - gamma Delta Vect(y)_(i,- 1), \
              Vect(Z)_i & = mat(
                            delim: "(", [ y_(i 0) ], 0, dots.h.c, 0 ;
                            0, [ y_(i 0),y_(i 1) ], dots.h.c, 0 ;
                            dots.v, dots.v, dots.down, dots.v;
                            0, 0, dots.h.c, [ y_(i 0),y_(i 1),...,y_(i,T - 2) ]
                          ),
$
то набор моментных условий вида @eq-ch7-AB.cond1 можно кратко записать как
$
  Ex (Vect(Z)^(prime)_i Delta Vect(epsilon)_i) = Ex [ Vect(Z)^(prime)_i (Delta Vect(y)_i - gamma Delta Vect(y)_(i,- 1)) ] = Vect(0)
$

Для получения оценки $gamma$ производится решение следующей задачи
$
  min_(gamma)[ 1/N sum^N_(i = 1)Vect(Z)^(prime)_i (Delta Vect(y)_i - gamma Delta Vect(y)_(i,- 1)) ]^(prime)Vect(W)_(N)[ 1/N sum^N_(i = 1)Vect(Z)^(prime)_i (Delta Vect(y)_i - gamma Delta Vect(y)_(i,- 1)) ] \
  hat(gamma)_(G M M) = [ (sum^N_(i = 1) Delta Vect(y)^(prime)_(i,- 1)Vect(Z)_i)^(prime)Vect(W)_N (sum^N_(i = 1)Vect(Z)^(prime)_i Delta Vect(y)_(i,- 1)) ]^(- 1)(sum^N_(i = 1) Delta Vect(y)^(prime)_(i,- 1)Vect(Z)_i)^(prime)Vect(W)_N (sum^N_(i = 1)Vect(Z)^(prime)_i Delta Vect(y)_i)
$

Если $Vect(W)_N$ положительно определена, то $hat(gamma)_(G M M)$ является состоятельной оценкой.
Для того, чтобы эта оценка была эффективной, матрица $Vect(W)_N$ должна быть асимптотически пропорциональной обратной ковариационной матрице выборочных моментных условий; оптимальной оценкой будет
$
  hat(Vect(W))^(o p t)_N = (1/N sum^N_(i = 1)Vect(Z)^(prime)_i Delta hat(epsilon)_i Delta hat(epsilon)^(prime)_i Vect(Z)_i)^(- 1),
$
где $hat(epsilon)_i$ — остатки модели с $Vect(W)_N = Vect(I)$.

Если постулировать гомоскедастичность и отсутствие автокорреляции, то
$
  Ex (Delta epsilon_i Delta epsilon_(i)') = sigma^2_(epsilon)Vect(G)& = sigma^2_(epsilon)mat(
    delim: "(", 2, - 1, 0, dots.h.c, 0 ;
    - 1, 2, - 1, dots.h.c, 0 ;
    0, - 1, 2, dots.h.c, 0 ;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    0, 0, 0, dots.h.c, 2
  ), \
  hat(Vect(W))^(o p t)_N & = (1/N sum^N_(i = 1)Vect(Z)^(prime)_i Vect(G) Vect(Z)_i)^(- 1),
$
то есть оценка может быть получена без предварительного этапа.

== Модель с экзогенными регрессорами

Модель с экзогенными переменными оценивается в целом аналогично модели без таковых.

Набор дополнительных инструментов зависит от природы добавляемых регрессоров.
Если регрессоры являются _строго экзогенными_, то есть $Ex (Vect(x)_(i s)epsilon_(k t)) = 0, forall i,j,s,t$, то
$ Ex (Vect(x)_(i s)Delta epsilon_(i t)) = 0, forall s,t ; $
это означает, что $Vect(x)_(i t)$ являются хорошими инструментами для _любого_ момента времени.

Однако добавление всех $Vect(x)_(i t)$ для всех $t$ приведет к очень большой матрице инструментов.
Поэтому в качестве инструментов рассматриваются первые разности экзогенных регрессоров, что дает следующие условия
$ Ex (Delta Vect(x)_(i t)Delta epsilon_(i t)) = 0, forall t $

Если регрессоры являются _предопределенными_, то есть некоррелированными с текущим и будущими значениями случайной ошибки, то $Vect(x)_(i,t - 1),...,Vect(x)_(i,1)$
будут подходящими инструментами для момента времени $t$ с моментными условиями
$ Ex (Vect(x)_(i,t - j)Delta epsilon_(i t)) = 0, j in { 1,...,t - 1 } $

Разумеется, что на практике часть экзогенных регрессоров будет строго экзогенными, а часть предопределенными, что потребует соответствующей
корректировки матриц $Vect(Z)_i$.

#bibliography("bibliography.bib")

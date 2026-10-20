#import "@preview/theorion:0.6.0": *
#import "@preview/tschich:0.2.0": *
#import "@local/neat-document:0.1.0": *
#import "@local/trinkets:0.1.0": *
#set document(title: [Классическая линейная модель и свойства МНК-оценок])
#let margins = tschich-var(210mm, 297mm, 1 / 10)
#let font-size = 14pt
#show: neat-document.with(margin: margins, leading: 1em, fonts: (size: font-size))
#show: show-theorion
#show: trinkets

#let TSS = math.italic("TSS")
#let ESS = math.italic("ESS")
#let RSS = math.italic("RSS")


= Модель парной линейной регрессии

Модель *парной линейной регрессии* имеет вид #cite(<nosko2021ekonometrika>)
#eqnum[$ y_i = beta_0 + beta_1 x_i + u_i $] <eq-2.1>
где
- $y_i$ — _зависимая_, или _объясняемая_, или _эндогенная_ переменная; очень редко встречается название _регрессанд_;
- $x_i$ — _объясняющая_, или _экзогенная_ переменная; очень часто используется название _регрессор_;
- $u_i$ — случайная ошибка.

Случайная ошибка включается в себя _все_ значимые, но не включенные в модель переменные.

Не стоит путать модель регрессии и *теоретическую линию регрессии*, имеющую вид
$ y= beta_0 + beta_1 x $
и показывающую теоретическую связь между переменными в рамках всей генеральной совокупности.

Оценить коэффициенты $beta_0$ и $beta_1$ проще всего при помощи метода наименьших квадратов, решая задачу
#eqnum[$ min_(beta_0,beta_1)underbrace(sum^N_(i = 1)(y_i - beta_0 - beta_1 x_i)^2, Q), $] <eq-2.2>
что приводит нас к условиям первого порядка
#eqnum[$
  frac(partial Q, partial beta_0)bar.v_(hat(beta)_0,hat(beta)_1) & = - 2 sum^N_(i = 1)(y_i - hat(beta)_0 - hat(beta)_1 x_i) = 0
$] <eq-2.3>
#eqnum[$
  frac(partial Q, partial beta_1)bar.v_(hat(beta)_0,hat(beta)_1) & = - 2 sum^N_(i = 1)(y_i - hat(beta)_0 - hat(beta)_1 x_i) x_i = 0
$] <eq-2.4>

Данные условия представляют собой систему из двух уравнений с двумя неизвестными, решение которого равно
$
  hat(beta)_1 & = frac(sum^N_(i = 1)(y_i - macron(y)) (x_i - macron(x)), sum^N_(i = 1)(x_i - macron(x))^2) = frac(s_(y x), s^2_x) \
  hat(beta)_0 & = macron(y) - hat(beta)_1 macron(x) \
  hat(u)_i & = y_i - hat(y)_i = y_i - hat(beta)_0 - hat(beta)_1 x_i
$

== Предпосылки МНК

+ Условное распределение $u_i$ относительно $X_i$ имеет нулевое среднее: $Ex [ u_i vert x_i ] = 0$;
  #image("fig/OLS_homo_errors.pdf", width: 8cm)
+ Наблюдения $(y_i,x_i) tilde #text[i.i.d.]$;
+ Большие выбросы маловероятны: $y_i$ и $x_i$ имеют ненулевые и конечные четвертые моменты.
  #image("fig/OLS_outliers.pdf", width: 8cm)

При выполнении этих предположений оценки коэффициентов $hat(beta)_0$ и $hat(beta)_1$ будут асимптотически нормальными.

== Графическая интерпретация

Для графической интерпретации МНК представим оценку @eq-2.1 в виде #cite(<nosko2021ekonometrika>)
$
  underbrace(
    mat(
      delim: "(", y_(1); y_(2);
      dots.v; y_N
    ), Vect(y)
  ) = hat(beta)_0 underbrace(
    mat(
      delim: "(", 1 ;
      1 ;
      dots.v;
      1
    ), Vect(1)
  ) + hat(beta)_1 underbrace(
    mat(
      delim: "(", x_(1); x_(2);
      dots.v; x_N
    ), Vect(x)
  ) + underbrace(
    mat(
      delim: "(", hat(u)_(1);
      hat(u)_(2);
      dots.v;
      hat(u)_N
    ), hat(Vect(u))
  ) = underbrace(hat(beta)_0 Vect(1) + hat(beta)_1 Vect(x), hat(Vect(y))) + hat(Vect(u)) = hat(Vect(y)) + hat(Vect(u))
$

Нетрудно увидеть, что условия первого порядка @eq-2.3 и @eq-2.4 примут вид
$
  lr(chevron.l hat(Vect(u)), Vect(1) chevron.r) & = 0, \
  lr(chevron.l hat(Vect(u)), Vect(x) chevron.r) & = 0,
$
то есть вектор остатков ортогонален векторам $upright(bold(1))$ и $x$; таким образом он ортогонален и "плоскости", натянутой на эти векторы.
Кроме того, во введенных обозначениях задача МНК представляет собой задачу поиска ортогональной проекции вектора $y$ на данную "плоскость".

Действительно, @eq-2.2 примет вид
$
  min_(beta_0,beta_1) autolr(Vect(y) underbrace(- hat(beta)_0 Vect(1)- hat(beta)_1 Vect(x), hat(Vect(y))), delim: #("\u{2016}", "\u{2016}"))^2,
$
из которой понятно, что $hat(y)$ — ортогональная проекция $y$ на "плоскость".
#align(center)[#image("fig/OLS.pdf")]

= Модель множественной линейной регрессии

Очевидно, что модель парной регрессии в большинстве случаев не может удовлетворить исследователя.
Более того, в такой модели крайне высока вероятность смещения из-за пропущенной переменной.

Модель *множественной линейной регрессии* имеет вид:
$ y_i = beta_0 + beta_1 x_(1 i) + beta_2 x_(2 i) + ...+ beta_K x_(K i) + u_i $

Оценки коэффициентов получаются аналогично парному случаю:
$ min_(beta_0,beta_1,...,beta_K)underbrace(sum^N_(i = 1)(y_i - beta_0 - beta_1 x_(1 i) - ...- beta_K x_(K i))^2, Q) $

Разумеется, что данную задачу можно решить "в лоб", найдя условия первого порядка в точке оптимума, приравняв их нулю и решив получившуюся систему уравнений.
Однако уже для случая $K = 2$ получающиеся выражения становятся весьма громоздкими.
Поэтому чаще всего работают с моделями в матричной форме #cite(<magnus2021ekonometrika>):
$
  underbrace(
    mat(
      delim: "(", y_(1); y_(2);
      dots.v; y_N
    ), Vect(y)
  )       & = underbrace(
              mat(
                delim: "(", 1, x_(11), x_(21), ..., x_(K 1);
                1, x_(12), x_(22), ..., x_(K 2);
                dots.v, dots.v, dots.v, dots.down, dots.v;
                1, x_(1 N), x_(2 N), ..., x_(italic("italic("K N")"))
              ), Vect(X)
            ) underbrace(
              mat(
                delim: "(", beta_(0);
                beta_(1);
                dots.v;
                beta_K
              ), Vect(beta)
            ) + underbrace(
              mat(
                delim: "(", u_(1); u_(2);
                dots.v; u_N
              ), Vect(u)
            ) \
  Vect(y) & = Vect(X) Vect(beta) + Vect(u)
$

Сумма квадратов остатков примет вид $hat(Vect(u))' hat(Vect(u))$, а оценка вектора коэффициентов будет равна
$ hat(Vect(beta)) = (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)Vect(y) $

#proof[
  Задача поиска оптимального вектора $hat(Vect(beta))$ имеет вид
  $ min_(Vect(beta))underbrace((Vect(Y) - Vect(X) Vect(beta))^(prime)(Vect(Y) - Vect(X) Vect(beta)), Vect(Q)) $

  Найдем дифференциал этого выражения #cite(<magnus2021ekonometrika>) #cite(<pierse2021ols>):
  $
    dif Vect(Q) & = dif [ (Vect(Y) - Vect(X) Vect(beta))^(prime) ] (Vect(Y) - Vect(X) Vect(beta)) + (Vect(Y) - Vect(X) Vect(beta))^(prime)dif (Vect(Y) - Vect(X) Vect(beta)) = \
    & = [ dif (Vect(Y) - Vect(X) Vect(beta))^(prime) ] (Vect(Y) - Vect(X) Vect(beta)) + (Vect(Y) - Vect(X) Vect(beta))^(prime)dif (Vect(Y) - Vect(X) Vect(beta)) = \
    & = (Vect(Y) - Vect(X) Vect(beta))^(prime)dif (Vect(Y) - Vect(X) Vect(beta)) + (Vect(Y) - Vect(X) Vect(beta))^(prime)dif (Vect(Y) - Vect(X) Vect(beta)) = \
    & = - 2(Vect(Y) - Vect(X) Vect(beta))^(prime)Vect(X) dif Vect(beta)
  $

  При выводе мы воспользовались тем, что $Vect(Q)$ — скаляр, а значит скалярами являются и слагаемые в промежуточных выражениях.
  А для скаляра, очевидно, справедливо $x^(prime) = x$.
  Также мы считаем матрицы $Vect(Y)$ и $Vect(X)$ фиксированными в рамках имеющейся выборки.

  Имеем, что в точке оптимума
  $ nabla Vect(Q) prop (Vect(Y) - Vect(X) Vect(beta))^(prime)Vect(X) = Vect(0) $
  откуда транспонированием и простыми преобразованиями и получаем искомую формулу.
]

Нетрудно увидеть, что подобранные значения равны $hat(Vect(y)) = Vect(X) (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)Vect(y) = Vect(H) Vect(y)$,
а остатки — $hat(Vect(u)) = Vect(y) - hat(Vect(y)) = [ Vect(I)_N - Vect(H) ] Vect(y) = Vect(M) Vect(y)$.
Также нетрудно увидеть, что
$
  Vect(H)^(prime) & = Vect(H) & Vect(M)^(prime) & = Vect(M) \
  Vect(H) Vect(H) & = Vect(H) & Vect(M) Vect(M) & = Vect(M) \
  Vect(M) Vect(H) & = Vect(0) & Vect(H) Vect(M) & = Vect(0) \
  Vect(H) Vect(X) & = Vect(X) & Vect(M) Vect(X) & = Vect(0)
$

Более того, так как матрица $Vect(M)$ идемпотентна, то $rank Vect(M) = tr Vect(M)$,
то есть
$
  tr Vect(M) & = tr (Vect(I)_N - Vect(X) (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)) = \
             & = tr Vect(I)_N - tr (Vect(X) (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)) = \
             & = N - underbrace(tr (Vect(X)^(prime)Vect(X) (Vect(X)^(prime)Vect(X))^(- 1)), Vect(I)_K) = N - K
$

Этот факт позволяет получить несмещенную оценку $sigma^2_u$ #cite(<pierse2021ols>):
$
  Ex (hat(Vect(u))^(prime)hat(Vect(u)) abs(thin Vect(X)) & = tr [ Ex (Vect(u)^(prime)Vect(M)^(prime)Vect(M) Vect(u)) thin Vect(X)) ] = Ex [ tr (Vect(u)^(prime)Vect(M) Vect(u)) vert thin Vect(X) ] = \
  & = Ex [ tr (Vect(M) Vect(u) Vect(u)^(prime)) vert thin Vect(X) ] = sigma^2_u tr Vect(M) = sigma^2_u (N - K) \
  sigma^2_u & = Ex (frac(hat(Vect(u))^(prime)hat(Vect(u)), N - K) vert thin Vect(X))
$

== Расширенные предпосылки МНК

+ $Ex (u_i divides Vect(X)_i) = 0$.
+ $(Vect(X)_i,Y_i) tilde #text[i.i.d.]$
+ $Vect(X)_i$ и $u_i$ имеют ненулевые конечные четвертые моменты.
+ $Vect(X)$ имеет полный ранг — нет совершенной мультиколлинеарности, то есть точной линейной зависимости между регрессорами.
+ $Var (u_i divides Vect(X)_i) = sigma^2_u$ — гомоскедастичность.
+ $u_i divides Vect(X)_i tilde.op N (0,sigma^2_u)$.

= Свойства МНК оценок

+ $hat(Vect(beta))$ — несмещенная оценка;
+ $display(frac(1, N))sum^N_(i = 1)hat(u)_i = 0$ — выборочное среднее остатков равно нулю;
+ $display(frac(1, N))sum^N_(i = 1)hat(y)_i = macron(y)$;
+ $sum^N_(i = 1)hat(u)_i x_(k i) = 0, forall k in { 1,...,K }$.
+ $TSS = ESS + RSS$, где
  - $TSS = sum^N_(i = 1)(y_i - macron(y))^2$ — полная сумма квадратов,
  - $ESS = sum^N_(i = 1)(hat(y_i) - macron(y))^2$ — объясненная сумма квадратов,
  - $RSS = sum^N_(i = 1)(y_i - hat(y)_i)^2 = sum^N_(i = 1)hat(u)^2_i$ — остаточная сумма квадратов или сумма квадратов остатков.
  #proof[
    $
      TSS & = sum^N_(i = 1)(y_i - macron(y))^2 = sum^N_(i = 1)[ (y_i - hat(y)_i) + (hat(y)_i - macron(y)) ]^2 = \
      & = sum^N_(i = 1)(y_i - hat(y)_i)^2 + sum^N_(i = 1)(hat(y)_i - macron(y))^2 + 2 sum^N_(i = 1)(y_i - hat(y)_i) (hat(y)_i - macron(y)) = \
      & = RSS + ESS + 2 sum^N_(i = 1)hat(u)_i hat(y)_i - 2 macron(y) sum^N_(i = 1)hat(u)_i
    $

    Последнее слагаемое равно нулю в силу второго свойства.
    Сумму в предпоследнем слагаемом можно представить в виде
    $ sum^N_(i = 1)hat(u)_i (hat(beta)_0 + hat(beta)_1 x_(1 i) + ...+ hat(beta)_K x_(K i)) $
    Раскрыв скобки и воспользовавшись вторым и четвертым свойствами, мы получим, что это слагаемое также равно нулю.
  ]

Важно отметить, что все перечисленные свойства, за исключением, разве что, четвертого, справедливы _только_ при наличии в модели свободного члена!
#align(center)[#image("fig/OLS_nocons.pdf")]

= Теорема Гаусса-Маркова

== Одномерный случай

#theorem[Гаусса-Маркова][
  Если для всех $i,j = 1,...,n$ выполняются условия Гаусса-Маркова:
  + $Ex (u_i vert X_1,...,X_n) = 0$,
  + $Var (u_i vert X_1,...,X_n) = sigma^2_u, 0 < sigma^2_(u)< infinity$,
  + $Ex (u_i u_j vert X_1,...,X_n) = 0, i != j$, то МНК оценка $hat(beta)_1$ является наилучшей (эффективной) линейной условно не смещенной оценкой (BLUE) #cite(<gujarati2009basic>).
]

- УГ-М (1) следует из предположений (1) и (2)
- УГ-М (2) следует из предположения (2) и предположения о гомоскедастичности ошибок
- УГ-М (3) следует из предположения (2)

#proof[
  Линейность:
  $
    sum^n_(i = 1)(X_i - macron(X)) & = 0 => \
    sum^n_(i = 1)(X_i - macron(X)) (Y_i - macron(Y)) & = sum^n_(i = 1)(X_i - macron(X)) Y_i => \
    hat(beta)_1 = frac(sum^n_(i = 1)(X_i - macron(X)) (Y_i - macron(Y)), sum^n_(i = 1)(X_i - macron(X))^2) & = frac(sum^n_(i = 1)(X_i - macron(X)) Y_i, sum^n_(i = 1)(X_i - macron(X))^2) = sum^n_(i = 1)hat(a)_i Y_i
  $

  Условная несмещенность:
  $
    Ex (hat(beta)_1 vert X_1,...) & = beta_1 + Ex (frac(sum^n_(i = 1)(X_i - macron(X)) u_i, sum^n_(i = 1)(X_i - macron(X))^2) vert thin X_1,...) = \
    & = beta_1 + frac(sum^n_(i = 1)(X_i - macron(X)) Ex (u_i vert X_1,...), sum^n_(i = 1)(X_i - macron(X))^2) = \
    & = beta_1
  $

  Эффективность: пусть $tilde(beta)_1$ — любая линейная условно несмещенная оценка $beta_1$, то есть $tilde(beta)_1 = sum^n_(i = 1)a_i Y_i$
  и $Ex tilde(beta)_1 = beta_1$.
  Тогда $sum^n_(i = 1)a_i = 0$, $sum^n_(i = 1)a_i X_i = 1$.
  $
    tilde(beta)_1 & = sum^n_(i = 1)a_i Y_i = beta_0 sum^n_(i = 1)a_i + beta_1 sum^n_(i = 1)a_i X_i + sum^n_(i = 1)a_i u_i => \
    Ex (tilde(beta)_1 vert X_1,...) & = beta_0 underbrace(sum^n_(i = 1)a_i, = 0) + beta_1 underbrace(sum^n_(i = 1)a_i X_i, = 1) = beta_1 \
    tilde(beta)_1 - beta_1 & = sum^n_(i = 1)a_i u_i => \
    Var (tilde(beta)_1 vert X_1,...) & = Ex [ (sum^n_(i = 1)a_i u_i)^2 vert thin X_1,...
    ] = \
    & = sum^n_(i = 1)sum^n_(j = 1)a_i a_j Ex (u_i,u_j vert X_1,...) = sigma^2_u sum^n_(i = 1)a^2_i
  $

  Пусть $a_i = hat(a)_i + d_i$.
  Тогда $sum^n_(i = 1)a^2_i = sum^n_(i = 1)hat(a)_i + 2 sum^n_(i = 1)hat(a)_i d_i + sum^n_(i = 1)d^2_i$.
  $
    sum^n_(i = 1)hat(a)_i d_i & = sum^n_(i = 1)frac((X_i - macron(X)) d_i, sum^n_(i = 1)(X_i - macron(X))^2) = \
    & = frac(sum^n_(i = 1)X_i d_i - macron(X) sum^n_(i = 1)d_i, sum^n_(i = 1)(X_i - macron(X))^2) = \
    & = frac((sum^n_(i = 1)X_i a_i - sum^n_(i = 1)X_i hat(a)_i) - macron(X) (sum^n_(i = 1)a_i - sum^n_(i = 1)hat(a)_i), sum^n_(i = 1)(X_i - macron(X))^2) = 0 \
    Var (tilde(beta)_1 vert X_1,...) & = sigma^2_u sum^n_(i = 1)a^2_i = \
    & = sigma^2_u sum^n_(i = 1)hat(a)^2_i + sigma^2_u sum^n_(i = 1)d^2_i = \
    & = Var (hat(beta)_1 vert X_1,...) + sigma^2_u sum^n_(i = 1)d^2_i => \
    Var (tilde(beta)_1 vert X_1,...) & - Var (hat(beta)_1 vert X_1,...) = sigma^2_u sum^n_(i = 1)d^2_(i)> 0
  $ #qedhere
]

== Многомерный случай

#theorem[Гаусса-Маркова][
  Пусть
  + $Ex (Vect(U) divides Vect(X)) = Vect(0)_N$.
  + $Ex (Vect(U) Vect(U)^(prime)divides Vect(X)) = sigma^2_u Vect(I)_N$.
  + $Vect(X)$ имеет полный столбцовый ранг.
  Тогда МНК оценка $hat(Vect(beta))$ является наилучшей (эффективной) линейной условно не смещенной оценкой (BLUE)#cite(<nosko2021ekonometrika>) #cite(<magnus2021ekonometrika>) #cite(<pierse2021ols>).
]
#proof[
  Линейность $hat(Vect(beta))$ очевидна.
  Действительно, так как
  $ hat(Vect(beta)) = (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)Vect(y) $
  то несложно увидеть, что каждый компонент $hat(Vect(beta))$ представляет собой линейную комбинацию $Y_i$ с весами, определяемыми матрицей $Vect(X)$.
  Также отметим, что
  $
    hat(Vect(beta)) = (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)(Vect(X) Vect(beta) + Vect(u)) = Vect(beta) + (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)Vect(u)
  $

  Несмещенность также несложно показать:
  $
    Ex (hat(Vect(beta)) divides Vect(X)) & = Ex [ (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)Vect(y) divides Vect(X) ] = \
    & = Ex [ (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)(Vect(X) Vect(beta) + Vect(u)) divides Vect(X) ] = \
    & = Vect(beta) + (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)underbrace(Ex (Vect(u) divides Vect(X)), = 0) = Vect(beta) quad #text[по 1-му условию]
  $

  Перейдем к эффективности.
  Пусть $tilde(Vect(beta))$ — *произвольная* линейная условно несмещенная оценка $Vect(beta)$.
  Ее можно представить в виде
  $ tilde(Vect(beta)) = (Vect(B) + Vect(C)) Vect(Y) $
  где $Vect(B) = (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)$.

  Из несмещенности $tilde(Vect(beta))$
  $
    Vect(beta) = Ex (tilde(Vect(beta)) divides Vect(X)) = (Vect(B) + Vect(C)) Ex (Vect(Y) divides Vect(X)) = (Vect(B) + Vect(C)) Vect(X) Vect(beta) = (Vect(I)_(K + 1) + Vect(C) Vect(X)) Vect(beta)
  $
  следует, что $Vect(C) Vect(X) = Vect(0)$, что означает, что $Vect(B) Vect(C)'= Vect(0)$.

  Также отметим, что
  $
    tilde(Vect(beta)) & = (Vect(B) + Vect(C)) Vect(Y) = (Vect(B) + Vect(C)) (Vect(X) Vect(beta) + Vect(u)) = \
    & = underbrace(Vect(B) Vect(X), = Vect(I)_(K + 1)) Vect(beta) + underbrace(Vect(C) Vect(X), = Vect(0)) Vect(beta) + (Vect(B) + Vect(C)) Vect(u) = \
    & = Vect(beta) + (Vect(B) + Vect(C)) Vect(u)
  $

  Ковариационная матрица $hat(Vect(beta))$ равна
  $
    Sigma_(hat(Vect(beta))) = Ex [ (hat(Vect(beta)) - Vect(beta)) (hat(Vect(beta)) - Vect(beta))^(prime)divides Vect(X) ] & = Ex [ (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)Vect(u) Vect(u)^(prime)Vect(X) (Vect(X)^(prime)Vect(X))^(- 1)divides Vect(X) ] = \
    & = (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)underbrace(Ex (Vect(u) Vect(u)^(prime)divides Vect(X)), Sigma_u = sigma^2_u Vect(I)_N) Vect(X) (Vect(X)^(prime)Vect(X))^(- 1) = \
    & = sigma^2_u (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)Vect(I)_N Vect(X) (Vect(X)^(prime)Vect(X))^(- 1) = \
    & = sigma^2_u (Vect(X)^(prime)Vect(X))^(- 1)
  $

  Ковариационная матрица $tilde(Vect(beta))$ равна
  $
    Sigma_(tilde(Vect(beta))) = Ex [ (tilde(Vect(beta)) - Vect(beta)) (tilde(Vect(beta)) - Vect(beta))^(prime)divides Vect(X) ] & = (Vect(B) + Vect(C)) underbrace(Ex (Vect(u) Vect(u)^(prime)divides Vect(X)), Sigma_u = sigma^2_u Vect(I)_N) (Vect(B)^(prime) + Vect(C)^(prime)) = \
    & = sigma^2_u (Vect(B) Vect(B)^(prime) + Vect(B) Vect(C)^(prime) + Vect(C) Vect(B)^(prime) + Vect(C) Vect(C)^(prime)) = \
    & = Sigma_(hat(Vect(beta))) + sigma^2_u Vect(C) Vect(C)^(prime)
  $

  Матрица $Vect(C) Vect(C)'$ — симметрична и положительно полуопределена, то есть $forall Vect(x)$
  $ Vect(x) Vect(C) Vect(C)^(prime)Vect(x) >= 0 $
  Это означает, что диагональные элементы матрицы неотрицательны.
  А так как на главных диагоналях ковариационных матриц располагаются	дисперсии оценок отдельных коэффициентов, то это означает, что
  $ forall i in {0,1,...,k }thick Var tilde(beta)_(i)>= Var hat(beta)_i $
]

= Критерии качества приближения данных моделью

/ Стандартная ошибка регрессии (SER):
$ italic("SER") = s_(hat(u)) $
где
$ s^2_(hat(u)) = frac(1, n-k-1)sum^n_(i = 1)hat(u)^2_i = frac(RSS, n-k-1) $
$italic("SER")$ — _standard error of the regression_ это оценка стандартного отклонения случайных ошибок регрессии $hat(u)_i$.

/ Коэффициент детерминации регрессии $R^2$:
$ R^2 = frac(ESS, TSS) = 1 - frac(RSS, TSS) = [ Corr (X,y) ]^2 $
$R^2$ — доля выборочной дисперсии $Y_i$, объясненной (или предсказанной) $X_i$.
$R^2$ всегда увеличивается при включении новых переменных в модель, так как $RSS$ всегда уменьшается.

/ Скорректированный коэффициент детерминации $macron(R)^2$:
$ macron(R)^2 = 1 - frac(n-1, n-k-1)frac(RSS, TSS) = 1 - frac(s^2_(hat(u)), s^2_Y) $

/ Информационные критерии:
$
  italic("AIC") & = 2 k - 2 ln hat(L) & italic("AIC") & = 2 k + N ln hat(sigma)^2 \
  italic("BIC") & = k ln N - 2 ln hat(L) & italic("BIC") & = underbrace(k ln N+ N ln hat(sigma)^2, #text[в линейной\ нормальной регрессии])
$

= Тестирование гипотез

== Тестирование единичной гипотезы

Рассмотрим модель
$ y_i = beta_0 + beta_1 x_(1 i) + beta_2 x_(2 i) + ...+ beta_K x_(K i) + u_i $
Предположим, что мы хотим протестировать одиночную гипотезу для какого-то коэффициента
$
  H_0 & : beta_j = beta_(j,0) \
  H_a & : beta_(j) != beta_(j,0)
$
где $j = 1,...,K$.
Для этого мы можем воспользоваться стандартным t-тестом.

+ Вычисляем стандартную ошибку $SE(hat(beta)_j))$.
+ Вычисляем тестовую статистику $t^"act" = frac(hat(beta)_j - beta_(j,0), SE(hat(beta)_j))$.
+ Отвергаем нулевую гипотезу на уровне значимости $alpha$, если $abs(t^"act") > t^"crit"_(1 - alpha / 2)$.
  Или, эквивалентно, если $p$-значение меньше $alpha$.
  Или, эквивалентно, если $hat(beta)_0$ не попадает в доверительный интервал ${ hat(beta)_j +- t^"crit"_(1 - alpha slash 2) times SE(hat(beta)_j) }$.

На больших выборках можно пользоваться критическими значениями для стандартного нормального распределения; на не очень больших выборках нужно использовать критические значения для распределения Стьюдента.

=== Тестирование более сложных гипотез

Допустим, что мы тестируем нулевую гипотезу $H_(0): Vect(w)'Vect(beta) = w_0$.

Известно, что
$
  Vect(w)^(prime)(hat(Vect(beta)) - Vect(beta)) & tilde.op N (0,sigma^2_u Vect(w)^(prime)(Vect(X)^(prime)Vect(X))^(-1) Vect(w)) \
  frac(Vect(w)^(prime)hat(Vect(beta)) - Vect(w)^(prime)Vect(beta), sigma_u sqrt(Vect(w)^(prime)(Vect(X)^(prime)Vect(X))^(- 1)Vect(w))) & tilde.op N(0,1)
$

Из формулы @eq-2.5 следует, что
$
  display(frac(frac(Vect(w)^(prime)hat(Vect(beta)) - Vect(w)^(prime)Vect(beta), sigma_u sqrt(Vect(w)^(prime)(Vect(X)^(prime)Vect(X))^(- 1)Vect(w))), sqrt(frac(1, N - K)(N - K) frac(hat(sigma)^2_u, sigma^2_u)))) = frac(Vect(w)^(prime)Vect(beta) - w_0, hat(sigma)_u sqrt(Vect(w)^(prime)(Vect(X)^(prime)Vect(X))^(- 1)Vect(w)))tilde.op t_(N - K)
$
Если $Vect(w) = Vect(e)_j = #autolr($0,...,0, underbrace(1, j),0,...,0$)'$,
то полученная формула примет вид
$ frac(hat(beta)_j, SE (hat(beta)_j))tilde.op t_(N - K) $

== Тестирование совместных гипотез

*Совместная гипотеза* накладывает условия (ограничения) на два (или более) коэффициента регрессии.
Первой мыслью может быть использование индивидуальных t-статистик для тестирования ограничений.
Рассмотрим вопрос, стоит ли отвергать ли совместную гипотезу на уровне 5%, если каждая индивидуальная $t$-статистика превышает 1.96 по
модулю.

Пусть $t_1$ и $t_2$ — индивидуальные $t$-статистики, для простоты независимо распределенные:
$
  t_1 = frac(hat(beta)_1 - 0, italic("SE")(hat(beta)_1)) #text[и] t_2 = frac(hat(beta)_2 - 0, italic("SE")(hat(beta)_2))
$

Предлагаемая выше процедура тестирования предлагает отвергать
$ H_(0): beta_1 = beta_2 = 0 $
если
$ abs(t_1) > 1.96 #text[и/или] abs(t_2) > 1.96 $

Вероятность некорректно отвергнуть $H_0$ на 5%-м уровне значимости, когда она верна, в такой процедуре равна
#fill-width-with-text[$
  Pr_(H_0)(abs(t_1) > 1.96 #text[и/или] abs(t_2) > 1.96) & = 1 - Pr_(H_0)(abs(t_1) <= 1.96 #text[и] abs(t_2) <= 1.96) = \
  & = 1 - Pr_(H_0)(abs(t_1) <= 1.96) times Pr_(H_0)(abs(t_2) <= 1.96) = \
  & = 1 - 0.95^2 = 0.0975 = 9.75 %
$]

- Предложенная процедура (на основе $t$-статистик) имеет неправильный размер критерия ($!= 5%$).
- Этот размер критерия зависит от корреляции между $t_1$ и $t_2$ и, следовательно, между $hat(beta)_1$ и $hat(beta)_2$.

Возможным решением является использование поправки Бонферрони: вместо $alpha$ в качестве уровня значимости следует использовать $display(frac(alpha, q))$.

=== Проверка совместных гипотез: F-статистика

$F$-статистика позволяет проверять одновременно все ограничения совместной гипотезы.
Для случая гомоскедастичных ошибок $F$-статистику можно вычислить следующим образом #cite(<nosko2021ekonometrika>):
- Оцените две регрессии: регрессию "с ограничениями" (т.е. в предположениях нулевой гипотезы) и регрессию "без ограничений" (в предположениях альтернативной).
- Вычислите $R^2_"R"$ и $R^2_"UR"$, соответственно, или $RSS_"R"$ и $RSS_"UR"$.
  Тогда
  $
    F = frac((R^2_"UR" - R^2_"R") slash q, (1-R^2_"UR") slash (n-k_"UR"-1)) = frac((italic("RSS")_"R" - italic("RSS")_"UR") slash q, italic("RSS")_"UR" slash (n-k_"UR"-1))
  $
  - $R^2_"R",R^2_"UR"$ — $R^2$ в регрессиях с ограничениями и без.
  - $italic("RSS")_"R",italic("RSS")_"UR"$ — RSS в регрессиях с ограничениями и без.
  - $q$ — число ограничений в нулевой гипотезе.
  - $k_"UR"$ — число объясняющих переменных (без константы) в модели без ограничений.

$F$-статистика для случая гомоскедастичных ошибок отвергает нулевую гипотезу, когда включение в модель двух (или больше) новых переменных увеличивает коэффициент детерминации "достаточно" для того, чтобы качество регрессии улучшилось "достаточно".
$F$-статистика для случая гомоскедастичных ошибок асимптотически распределена как $chi^2_q slash q$#footnote[Если ошибки на самом деле гетероскедастичны, то асимптотическое распределение $F$-статистики для случая гомоскедастичных ошибок не является $chi^2_q slash q$.]:
#align(center)[
  #table(
    columns: (auto, auto),
    stroke: none,
    align: (center, center),
    table.hline(),
    [$q$], [5%],
    table.hline(),
    [1], [3.84],
    [2], [3.00],
    [3], [2.60],
    [4], [2.37],
    [5], [2.21],
    table.hline(),
  )
]

=== F-статистика в матричной форме

Допустим, что мы тестируем нулевую гипотезу $H_(0): Vect(R) Vect(beta) = Vect(r)$.

Известно, что
$
  Vect(R) hat(Vect(beta)) - Vect(r) = Vect(R) (hat(Vect(beta)) - Vect(beta)) tilde.op N (0,sigma^2_u Vect(R) (Vect(X)^(prime)Vect(X))^(- 1)Vect(R)^(prime))
$

Существует матрица $Vect(Q)$ такая, что $[ Vect(R) (Vect(X)'Vect(X))^(- 1)Vect(R)' ]^(- 1) = Vect(Q)'Vect(Q)$, и что
$ Vect(xi) = Vect(Q) Vect(R) (hat(Vect(beta)) - Vect(beta)) tilde.op N (0,Vect(I)_N) $

Поэтому при справедливости нулевой гипотезы #cite(<magnus2021ekonometrika>) #cite(<pierse2021ols>)
$
  frac(Vect(xi)' Vect(xi), sigma^2_u) & = frac((Vect(R) hat(Vect(beta)) - Vect(r))^(prime)Vect(Q)^(prime)Vect(Q) Vect(R) (Vect(R) hat(Vect(beta)) - Vect(r)), sigma^2_u) = \
  & = frac((Vect(R) hat(Vect(beta)) - Vect(r))^(prime)[ Vect(R) (Vect(X)^(prime)Vect(X))^(- 1)Vect(R)^(prime) ]^(- 1)(Vect(R) hat(Vect(beta)) - Vect(r)), sigma^2_u) = \
  & = frac((hat(Vect(beta)) - Vect(beta))^(prime)Vect(R)^(prime)[ Vect(R) (Vect(X)^(prime)Vect(X))^(- 1)Vect(R)^(prime) ]^(- 1)Vect(R) (hat(Vect(beta)) - Vect(beta)), sigma^2_u) = \
  & = frac(Vect(u)^(prime)Vect(X) (Vect(X)^(prime)Vect(X))^(- 1)Vect(R)^(prime)[ Vect(R) (Vect(X)^(prime)Vect(X))^(- 1)Vect(R)^(prime) ]^(- 1)Vect(R) (Vect(X)^(prime)Vect(X))^(- 1)Vect(X)^(prime)Vect(u), sigma^2_u)tilde chi^2_q
$
Это следует из идемпотентности матрицы в числителе.

С другой стороны,
$
  frac(hat(Vect(u))^(prime)hat(Vect(u)), sigma^2_u) = frac(Vect(u)^(prime)Vect(M) Vect(u), sigma^2_u) = frac((N - K) hat(sigma)^2_u, sigma^2_u)tilde chi^2_(N - K)
$ <eq-2.5>
откуда получим
$
  F = display(frac(1/q frac(Vect(xi)' Vect(xi), sigma^2_u), frac(1, N - K)frac((N - K) hat(sigma)^2_u, sigma^2_u))) = display(frac((Vect(R) hat(Vect(beta)) - Vect(r))^(prime)[ Vect(R) (Vect(X)^(prime)Vect(X))^(- 1)Vect(R)^(prime) ]^(- 1)(Vect(R) hat(Vect(beta)) - Vect(r)), q hat(sigma)^2_u))
$

Если ошибки гетероскедастичны, то асимптотическое распределение $F$-статистики для случая гомоскедастичных ошибок не является $chi^2_q slash q$.
Устойчивая к гетероскедастичности $F$-статистика (статистика Вальда)
имеет вид
$
  F & = frac((Vect(R) hat(Vect(beta)) - Vect(r))^(prime)(Vect(R) hat(Vect(Sigma))_(hat(beta))Vect(R)^(prime))^(- 1)(Vect(R) hat(Vect(beta)) - Vect(r)), q)tilde chi^2_q slash q \
  W & = (Vect(R) hat(Vect(beta)) - Vect(r))^(prime)(Vect(R) hat(Vect(Sigma))_(hat(beta))Vect(R)^(prime))^(- 1)(Vect(R) hat(Vect(beta)) - Vect(r)) = q F tilde chi^2_q
$
где $Vect(R)$ — матрица ограничений, $Vect(r)$ — вектор свободных членов ограничений, $hat(Vect(Sigma))_(hat(beta))$ — ковариационная матрица вектора $hat(Vect(beta))$.

Для того, чтобы $F$-статистика была _точной_, нужно, чтобы:
+ ошибки $u_i$ были гомоскедастичны, т.е. $Var (u_i vert X_1,...,X_k) = "const"$.
+ ошибки $u_1,...,u_n$ были нормально распределены.
В случае ненормальности ошибок регрессии статистика будет _асимптотической_.

Исследователи часто используют $F$-распределение на небольших выборках, так как оно имеет большие критические значения и является более консервативным.

= Фиктивные регрессоры

Пусть $X_i$ обозначает бинарную переменную.
Рассмотрим регрессию $Y_i = beta_0 + beta_1 X_i + u_i$.
Пусть $macron(Y)_0$ означает выборочное среднее для наблюдений с $X = 0$, а $macron(Y)_1$ обозначает
выборочное среднее для наблюдений с $X = 1$.
Покажем, что $hat(beta)_0 = macron(Y)_0$, $hat(beta)_0 + hat(beta)_1 = macron(Y)_1$ и $hat(beta)_1 = macron(Y)_1 - macron(Y)_0$.
$
                                  sum X_i & = n_1 \
                                macron(X) & = frac(n_1, n) \
                  sum (X_i - macron(X))^2 & = sum X^2_i - n macron(X)^2 = n_1 - n frac(n^2_1, n^2) - frac(n_0 n_1, n) \
  sum (X_i - macron(X)) (Y_i - macron(Y)) & = sum X_i (Y_i - macron(Y)) = sum X_i Y_i - macron(Y) sum X_i = \
                                          & = n_1 macron(Y)_1 - n_1 macron(Y)
$
$
  hat(beta)_1 & = frac(n_1 (macron(Y)_1 - macron(Y)), frac(n_0 n_1, n)) = frac(n, n_0)(macron(Y)_1 - macron(Y)) = frac(n, n_0)(macron(Y)_1 - frac(n_1, n)macron(Y)_1 - frac(n_0, n)macron(Y)_0) = \
  & = frac(n, n_0)(frac(n_0, n)macron(Y)_1 - frac(n_0, n)macron(Y)_0) = macron(Y)_1 - macron(Y)_0 \
  hat(beta)_0 & = macron(Y) - hat(beta)_1 macron(X) = macron(Y) - (macron(Y)_1 - macron(Y)_0) frac(n_1, n) = \
  & = frac(n_1, n)macron(Y)_1 + frac(n_0, n)macron(Y)_0 - frac(n_1, n)macron(Y)_1 - frac(n_1, n)macron(Y)_0 = macron(Y)_0
$

= Нелинейные модели

$ Y_i = f(X_(1 i),...,X_(k i)) + u_i, i = 1,...,n $
$f(X_(1 i),...,X_(k i))$ — *теоретическая функция нелинейной регрессии*.
Например:
$ f(italic("Income")_i) = beta_0 + beta_1 italic("Income")_i + beta_2 italic("Income") $

== Величина изменения $Y$

$ Delta Y = f(X_1 + Delta X_1,X_2,...,X_k) - f(X_1,X_2,...,X_k) $
Пусть $hat(f) (X_1,X_2,...,X_k)$ — предсказанное значение
$Y$, полученное на основе оценки $hat(f)$ функции теоретической
регрессии.
Тогда
$ Delta hat(Y) = hat(f) (X_1 + Delta X_1,X_2,...,X_k) - hat(f) (X_1,X_2,...,X_k) $

== Логарифмические регрессии

Дальнейшие выкладки опираются на следующее:
$ ln (x + Delta x) - ln (x)approx.eq frac(Delta x, x), #text[если] frac(Delta x, x) #text[мало] $
$frac(Delta x, x)$ — процентное изменение $x$, деленное на 100.

/ Линейно-логарифмическая модель: 1%-ное изменение $X$ приводит к изменению $Y$ на $display(frac(1, 100))beta_1$.
  $ Y_i = beta_0 + beta_1 ln (X_i) + u_i, i = 1,...,n $

/ Логарифмически-линейная модель: Единичное изменение $X$ (то есть $Delta X = 1$) приводит к изменению $Y$ на $100 beta_(1)%$.
  $ ln Y_i = beta_0 + beta_1 X_i + u_i, i = 1,...,n $
  Важно понимать, что так как $Delta ln Y= ln display(frac(Y + Delta Y, Y))$ лишь примерно равно $display(frac(Delta Y, Y))$, то приведенная выше интерпретация справедлива только в случае _малых_ значений $beta_1$.

  #table(
    columns: (auto,) * 4,
    stroke: none,
    align: (center,) * 4,
    table.hline(),
    [$Delta ln Y$],
    [$exp { Delta ln Y } = frac(Y + Delta Y, Y)$],
    [$Delta ln Y$],
    [$exp { Delta ln Y } = frac(Y + Delta Y, Y)$],
    table.hline(),
    [0.01], [1.01005], [0.1], [1.10520],
    [0.02], [1.02020], [0.2], [1.22140],
    [0.03], [1.03046], [0.3], [1.35000],
    [0.05], [1.05127], [0.5], [1.65000],
    [0.09], [1.09417], [1.0], [2.71828],
    table.hline(),
  )

/ Линейная в логарифмах: 1%-ное изменение $X$ (то есть $Delta X= 1$) приводит к изменению $Y$ на $beta_(1)%$.
$ ln Y_i = beta_0 + beta_1 ln X_i + u_i, i = 1,...,n $

#align(center)[
  #table(
    columns: (auto, auto, auto),
    stroke: none,
    align: (center, center, center),
    table.hline(),
    [Модель регрессии\ по $X$], [Теоретическая\ модель], [Эластичность\ $Ex (Y vert X)$],
    table.hline(),
    [Линейная], [$Y = beta_0 + beta_1 X + u$], [$frac(beta_1 Delta X, beta_0 + beta_1 X)$],
    [Линейно-логарифмическая], [$Y = beta_0 + beta_1 ln X+ u$], [$frac(beta_1, beta_0 + beta_1 ln X)$],
    [Логарифмически-линейная], [$ln Y= beta_0 + beta_1 X + u$], [$beta_1 Delta X$],
    [Линейная в логарифмах], [$ln Y= beta_0 + beta_1 ln X+ u$], [$beta_1$],
    table.hline(),
  )
]

*Используя критерии качества подгонки моделей, можно сравнивать модели, у которых одинаковый $Y$!*
- Можно сравнивать логарифмически-линейную и линейную в логарифмах модель.
- Можно сравнивать линейно-логарифмическую и полиномиальные модели.

#bibliography("bibliography.bib")

#show: appendix

= Регрессия по частям

Этот дополнительный материал будет полезным для рассмотрения некоторых последующих тем.
Смысл данного метода в том, что мы оцениваем коэффициенты модели не разом, а по очереди.

Рассмотрим модель #cite(<pierse2021gls>)
$
  Vect(y) = Vect(X) Vect(beta) + Vect(u) = Vect(X)_1 Vect(beta)_1 + Vect(X)_2 Vect(beta)_2 + Vect(u) = mat(delim: "(", Vect(X)_1, Vect(X)_2) mat(
    delim: "(", Vect(beta)_(1);
    Vect(beta)_2
  ) + Vect(u)
$
где $Vect(X)_1$ и $Vect(X)_2$ представляют собой разбиение матрицы
$Vect(X)$ по столбцам, причем столбцы в этих матрицах (а) не пересекаются, и (б) включают в себя все столбцы исходной матрицы.
Для дальнейшего повествования нам понадобятся формулы для матрицы, обратной к блочной:
#eqnum[$
  mat(delim: "(", Vect(A), Vect(B); Vect(C), Vect(D))^(- 1) & = mat(
                                                                delim: "(", tilde(Vect(A)), - tilde(Vect(A)) Vect(B) Vect(D)^(- 1);
                                                                - Vect(D)^(- 1)Vect(C) tilde(Vect(A)), Vect(D)^(- 1) + Vect(D)^(- 1)Vect(C) tilde(Vect(A)) Vect(B) Vect(D)^(- 1)
                                                              ) =
$] <eq-2.6>
#eqnum[$
  & = mat(delim: "(", Vect(A)^(- 1) + Vect(A)^(- 1)Vect(B) tilde(Vect(D)) Vect(C) Vect(A)^(- 1), - Vect(A)^(- 1)Vect(B) tilde(Vect(D)); - tilde(Vect(D)) Vect(C) Vect(A)^(- 1), tilde(Vect(D)))
$] <eq-2.7>
где $tilde(Vect(A)) = (Vect(A) - Vect(B) Vect(D)^(- 1)Vect(C))^(- 1)$, $tilde(Vect(D)) = (Vect(D) - Vect(C) Vect(A)^(- 1)Vect(B))^(- 1)$.

Очевидно, что формула @eq-2.6 справедлива в случае, когда матрица $Vect(D)$ обратима, а формула @eq-2.7 — когда обратима матрица $Vect(A)$.

Оценка $hat(Vect(beta))$ равна
$
  mat(
    delim: "(", hat(Vect(beta))_(1);
    hat(Vect(beta))_2
  ) = (Vect(X)'Vect(X))^(- 1)mat(
    delim: "(", Vect(X)_1 Vect(y);
    Vect(X)_2 Vect(y)
  )
$

Найдем $hat(Vect(beta))_1$.
Так как
$
  Vect(X)'Vect(X) = mat(
    delim: "(", Vect(X)^(prime)_1 Vect(X)_1, Vect(X)^(prime)_1 Vect(X)_(2);
    Vect(X)^(prime)_2 Vect(X)_1, Vect(X)^(prime)_2 Vect(X)_2
  )
$
то
$
  (Vect(X)^(prime)Vect(X))^(- 1) = mat(
    delim: "(", tilde(Vect(X)), - tilde(Vect(X)) Vect(X)^(prime)_1 Vect(X)_2 (Vect(X)^(prime)_2 Vect(X)_2)^(- 1);
    - (Vect(X)^(prime)_2 Vect(X)_2) Vect(X)^(prime)_2 Vect(X)_1 tilde(Vect(X)), (Vect(X)^(prime)_2 Vect(X)_2)^(- 1) + (Vect(X)^(prime)_2 Vect(X)_2)^(- 1)Vect(X)^(prime)_2 Vect(X)_1 tilde(Vect(X)) Vect(X)^(prime)_1 Vect(X)_2 (Vect(X)^(prime)_2 Vect(X)_2)^(- 1)
  )
$
где
$
  tilde(Vect(X)) & = [ Vect(X)^(prime)_1 Vect(X)_1 - Vect(X)^(prime)_1 Vect(X)_2 (Vect(X)^(prime)_2 Vect(X)_2)^(- 1)Vect(X)^(prime)_2 Vect(X)_1 ]^(- 1) = \
  & = { Vect(X)^(prime)_1 underbrace([ Vect(I)_N - Vect(X)_2 (Vect(X)^(prime)_2 Vect(X)_2)^(- 1)Vect(X)^(prime)_2 ], Vect(M)_2) Vect(X)^(prime)_1 }^(- 1)
$

Таким образом
$
  hat(Vect(beta))_1 & = tilde(Vect(X)) Vect(X)^(prime)_1 Vect(y) - tilde(Vect(X)) Vect(X)^(prime)_1 Vect(X)_2 (Vect(X)^(prime)_2 Vect(X)_2)^(- 1)Vect(X)^(prime)_2 Vect(y) = \
  & = tilde(Vect(X)) Vect(X)^(prime)_(1)[ Vect(I)_N - Vect(X)_2 (Vect(X)^(prime)_2 Vect(X)_2)^(- 1)Vect(X)^(prime)_2 ] Vect(y) = \
  & = (Vect(X)^(prime)_1 Vect(M)_2 Vect(X)_1)^(- 1)Vect(X)^(prime)_1 Vect(M)_2 Vect(y) = \
  & = (Vect(X)^(prime)_1 Vect(M)_2 Vect(M)_2 Vect(X)_1)^(- 1)Vect(X)^(prime)_1 Vect(M)_2 Vect(M)_2 Vect(y) = \
  & = (hat(Vect(X))^(prime)_1 hat(Vect(X))_1)^(- 1)hat(Vect(X))^(prime)_1 hat(Vect(y))
$
где $hat(Vect(X))_1$ и $hat(Vect(y))$ — остатки регрессий $Vect(X)_1$ и $Vect(y)$ на $Vect(X)_2$, соответственно.

Аналогично показывается, что
$ hat(Vect(beta))_2 = (Vect(X)^(prime)_2 Vect(M)_1 Vect(X)_2)^(- 1)Vect(X)^(prime)_2 Vect(M)_1 Vect(y) $

= Распределение хи-квадрат

Рассмотрим вектор $Vect(x)$ длины $p$ такой, что $Vect(x) tilde.op N (Vect(0),Vect(I)_p)$.
Тогда случайная величина $Vect(x)'Vect(x)$ будет иметь распределение хи-квадрат с $p$ степенями свободы $chi^2_p$.

Более того, если $Vect(x) tilde.op N (Vect(0),q^2 Vect(I)_p)$, то
$ frac(Vect(x)^(prime)Vect(x), q^2)tilde chi^2_p $

Пусть $Vect(A)$ — симметричная идемпотентная матрица ранга $r$.
Тогда #cite(<pierse2021gls>)
$ frac(Vect(x)^(prime)Vect(A) Vect(x), q^2)tilde chi^2_r $

#proof[
  Матрица $Vect(A)$ может быть спектрально разложена:
  $ Vect(A) = Vect(V) Vect(Lambda) Vect(V)^(prime) $
  где $Vect(V)$ — матрица собственных векторов такая, что $Vect(V)^(prime)Vect(V) = Vect(I)_p$, $Vect(Lambda)$ — диагональная матрица с собственными значениями на главной диагонали.

  Отсюда
  $
    Vect(x)^(prime)Vect(A) Vect(x) = Vect(x)^(prime)Vect(V) Vect(Lambda) Vect(V)^(prime)Vect(x) = Vect(z)^(prime)Vect(Lambda) Vect(z)
  $
  причем $Vect(z) tilde.op N (Vect(0),q^2 Vect(V)^(prime)Vect(V)) = N(Vect(0),q^2 Vect(I)_p)$.
  С другой стороны, так как $Vect(A)$ симметрична и идемпотентна, то ее собственные значения могут принимать значения только $0$ или $1$, причем число значений, равных $1$, равно рангу матрицы.

  Таким образом
  $
    Vect(z)^(prime)Vect(Lambda) Vect(z) = sum^p_(i = 1)lambda_i z^2_i = sum^r_(i = 1)z^2_i = Vect(z)^(prime)_1 Vect(z)_1
  $
  откуда следует
  $ frac(Vect(x)^(prime)Vect(A) Vect(x), q^2) = frac(Vect(z)^(prime)_1 Vect(z)_1, q^2)tilde chi^2_r $
]

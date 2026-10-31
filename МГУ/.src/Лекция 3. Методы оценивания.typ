#import "@preview/theorion:0.6.0": *
#import "@preview/tschich:0.2.0": *
#import "@local/neat-document:0.1.0": *
#import "@local/trinkets:0.1.0": *
#set document(title: [Методы оценивания], author: "Зямалов В.Е.")
#let margins = tschich-var(210mm, 297mm, 1 / 10)
#let font-size = 14pt
#show: neat-document.with(margin: margins, leading: 1em, fonts: (size: font-size))
#show: show-theorion
#show: trinkets

Основные источники: учебники Хаяши #cite(<hayyashi2017ekonometrika>) #cite(<hayashi2000econometrics>) и Гамильтона #cite(<hamilton1994time>).

Хаяши выделяет два класса так называемых "экстремальных" оценок:
- M-оценки, получаемые из задачи, целевая функция которой имеет вид
  $ Q_T (theta) = 1/T sum^T_(t = 1) m(y_t vert theta) $
  где $m(y_t vert theta)$ — некоторая функция $Vect(Y) times Theta -> RR$.
- GMM, получаемые из задачи, целевая функция которой имеет вид
  $ Q_T (theta) = - 1/2Vect(g)_T (theta)'hat(Vect(W)) Vect(g)_T (theta) $
  где $Vect(g)_T (theta) = 1/T sum^T_(t = 1) Vect(g) (y_t vert theta)$.

= Метод максимального правдоподобия

Рассмотрим следующую задачу.
У нас есть некоторый набор наблюдений ${y_t}$, полученный из генеральной совокупности с набором некоторых параметров $theta$, которые нам нужно оценить.

Для решения этой задачи можно воспользоваться *методом максимального правдоподобия*.
Предположим, что наблюдения являются независимыми и одинаково распределенными, причем вид распределения нам *известен*#footnote[Это утверждение приведено в учебнике #cite(<hayyashi2017ekonometrika>), хотя, как мы увидим далее, это не обязательное требование.
  Главное, чтобы мы *знали* распределения отдельных наблюдений.].
Для данных наблюдений можно записать плотность распределения каждого
из них с учетом значений вектора параметров $theta$ — $f(y_(t);theta)$.

Так как все наблюдения считаются независимыми, то достаточно легко можно вычислить совместную плотность распределения всего набора наблюдений
$ f(y_1, y_2, ..., y_T vert theta) = product^T_(t = 1) f(y_t vert theta) $

Полученное выражение и задает нам *функцию правдоподобия*.
Формально можно записать, что
$ L(theta vert Vect(y)) = Pr(Vect(y) vert theta) $
что означает, что значение функции правдоподобия, вычисленной на данном наборе параметров с учетом текущей выборки, численно совпадает с вероятностью получить данную выборку при данном наборе параметров.

Оптимальной оценкой $theta$
$ arg max_(theta in Theta) L(theta vert Vect(y)) $

Для упрощения функцию правдоподобия часто берут в логарифмах
$ l(theta vert Vect(y)) = sum^T_(t = 1) ln f(y_t vert theta) $
Необходимым условием экстремума для решения задачи является равенство нулю градиента функции правдоподобия
$
  nabla l(theta vert Vect(y)) bar.v_(theta = hat(theta)) = frac(italic("tiall"), partial theta) (hat(theta)) = Vect(0)
$

Достаточным условием является отрицательная определенность матрицы Гессе
$
  H(theta vert Vect(y)) = (frac(partial^2 l, partial theta partial theta')) (theta vert Vect(y))
$

Как отмечает Хаяши #cite(<hayyashi2017ekonometrika>), при выполнении достаточно простых условий ML-оценки являются состоятельными:

- функция правдоподобия должна быть вогнутой на пространстве параметров $Theta$;
- функция правдоподобия измерима на пространстве данных для любого $theta in Theta$.
- стационарный процесс ${y_t}$ эргодичен.

== Линейная регрессия

Рассмотрим пример оценивания линейной регрессии #cite(<magnus2021ekonometrika>).

Допустим, что мы ищем оценку параметров модели
$ y_t = beta_0 + beta_1 x_t + epsilon_t $

Предположим, что выполняются все предпосылки МНК, а также дополнительное предположение гауссовости случайных ошибок
$ epsilon_t sim "iid" N(0, sigma^2) $
что эквивалентно#footnote[В примере подразумевается детерминированность $x_t$.]
$ y_t sim N(beta_0 + beta_1 x_t, sigma^2) $

Для решения данной задачи имеет смысл рассмотреть так называемый *условный метод максимального правдоподобия*. Смысл метода заключается в том, что исходный набор наблюдений разбивается на две части, причем исследователя  интересует влияние одной части данных — в нашем случае $x_t$ — на условное матожидание других — $y_t$.
Плотность можно разложить следующим образом
$ f(y_t, x_t vert theta, psi) = f(y_t vert x_(t);theta) f(x_(t);psi) $
где $psi$ — параметры распределения $x_t$.

В итоге функция правдоподобия примет вид
$ l(theta vert y_t, x_t) = 1/T sum^T_(t = 1) ln f(y_t vert x_(t); theta) + 1/T sum^T_(t = 1) f(x_(t);psi) $
Очевидно, что при поиске оптимального значения $hat(theta)$ максимизироваться будет только первая часть; вторая часть игнорируется.

Вернемся к примеру.
Запишем функцию правдоподобия
$
  L(beta_0,beta_1,sigma^2 vert y_1,...,y_T) = Pr(y_1, ..., y_T vert x_1, ..., x_T, beta_0, beta_1, sigma^2) = \
  = (2 pi)^(- T slash 2) (sigma^2)^(- T slash 2) exp (- frac(1, 2 sigma^2) sum^T_(t = 1) (y_t - beta_0 - beta_1 x_t)^2)
$
или
$
  & ln L(beta_0, beta_1, sigma^2 vert y_1, ..., y_T) = \
  & = - T/2ln (2 pi) - T/2ln (sigma^2) - frac(1, 2 sigma^2) sum^T_(t = 1) (y_t - beta_0 - beta_1 x_t)^2
$

Найдем необходимые условия первого порядка максимизации функции
$
  frac(partial ln L, partial beta_0) & = frac(1, sigma^2) sum^T_(t = 1) (y_t - beta_0 - beta_1 x_t) = 0 \
  frac(partial ln L, partial beta_1) & = frac(1, sigma^2) sum^T_(t = 1) x_t (y_t - beta_0 - beta_1 x_t) = 0 \
  frac(partial ln L, partial sigma^2) & = - n/2frac(1, sigma^2) + frac(1, 2 sigma^4) sum^T_(t = 1) (y_t - beta_0 - beta_1 x_t)^2 = 0
$
Нетрудно заметить, что данные уравнения дают нам оценки
$
   hat(beta)_0 & = overline(y) - hat(beta)_1 overline(x) \
   hat(beta)_1 & = frac(sum (x_t - overline(x)) (y_t - overline(y)), sum (x_t - overline(x))^2) \
  hat(sigma)^2 & = 1/n sum e^2_t
$
То есть полученные оценки коэффициентов совпадают с МНК оценками.
При этом оценка $sigma^2$ оказывается состоятельной, но смещенной.

== AR(p)

Рассмотрим функцию правдоподобия для процесса AR($p$)
$ Y_t = c + alpha_1 Y_(t - 1) + ...+ alpha_p Y_(t - p) + epsilon_t $
где $epsilon sim "WN"(0, sigma^2)$.
Оцениваемый вектор параметров
$theta = (c,alpha_1,...,alpha_p,sigma^2)$.

Очевидно, что с описанием данного процесса есть сложность, связанная с первыми $p$ наблюдениями.
Гамильтон #cite(<hamilton1994time>) использует способ разделения функции на две части.
Первая часть определяет первые $p$ наблюдений как реализацию $p$-мерной Гауссовой случайной величины с математическим ожиданием
$ mu = frac(c, 1 - sum^p_(i = 1) alpha_i) $
Ковариационная матрица $(y_1,...,y_p)$ равна
$
  sigma^2 Vect(V)_p & = mat(
                        delim: "[", Ex (Y_1 - mu)^2, Ex (Y_1 - mu) (Y_2 - mu), dots.h.c, Ex (Y_1 - mu) (Y_p - mu);
                        Ex (Y_2 - mu) (Y_1 - mu), Ex (Y_2 - mu)^2, dots.h.c, Ex (Y_2 - mu) (Y_p - mu);
                        dots.v, dots.v, dots.down, dots.v;
                        Ex (Y_p - mu) (Y_1 - mu), Ex (Y_p - mu) (Y_2 - mu), dots.h.c, Ex (Y_p - mu)^2
                      ) = \
                    & = mat(
                        delim: "[", gamma_0, gamma_1, dots.h.c, gamma_(p - 1);
                        gamma_1, gamma_0, dots.h.c, gamma_(p - 2);
                        dots.v, dots.v, dots.down, dots.v;
                        gamma_(p - 1), gamma_(p - 2), dots.h.c, gamma_0
                      )
$
Отсюда совместная плотность первых $p$ наблюдений равна
$
  & f_(Y_p,Y_(p - 1),...,Y_1) (y_p,...,y_(1);theta) = \
  & wide = (2 pi)^(- p slash 2) | sigma^(- 2) Vect(V)^(- 1)_(p) |^(1 slash 2) exp [ - frac(1, 2 sigma^2) (Vect(y)_p - Vect(mu)_p)'Vect(V)^(- 1)_p (Vect(y)_p - Vect(mu)_p) ] = \
  & wide = (2 pi)^(- p slash 2) (sigma^(- 2))^(p slash 2) | Vect(V)^(- 1)_(p) |^(1 slash 2) exp [ - frac(1, 2 sigma^2) (Vect(y)_p - Vect(mu)_p)'Vect(V)^(- 1)_p (Vect(y)_p - Vect(mu)_p) ]
$
#cite(<galbraith1974inverses>) показали, что элемент $i j$ матрицы $Vect(V)^(- 1)_p$ равен
$
  v^(i j) = v^(j i) = [ sum^(i - 1)_(k = 0) alpha_k alpha_(k + j - i) - sum^(p + i - j)_(k = p + 1 - j) alpha_k alpha_(k + j - i) ]
$
причем $alpha_0 = - 1$#footnote[В исходном тексте свободный член в модели обозначается как $c$, чтобы избежать коллизии с $v^(i j)$.].

Для оставшейся части выборки поступают следующим образом.
Для наблюдения $t$ математическое ожидание зависит *только* от предыдущих $p$ наблюдений, то есть
$
  Ex (Y_t vert Y_(t - 1),...,Y_1) = Ex (Y_t vert Y_(t - 1),...,Y_(t - p)) = alpha_0 + sum^p_(i = 1) alpha_i y_(t - i)
$
При этом дисперсия равна $sigma^2$.
Поэтому для второй части плотность $Y_t$ имеет вид
$
  & f_(Y_t vert Y_(t - 1),...,Y_(t - p)) (y_t vert y_(t - 1),...,y_(t - p);theta) = \
  & wide = frac(1, sqrt(2 pi sigma^2)) exp [ - frac((y_t - alpha_0 - sum^p_(i = 1) alpha_i y_(t - i))^2, 2 sigma^2) ]
$
А совокупная плотность *всей* выборки равна
$
  & f_(Y_T,...,Y_1) (y_T,...,y_(1);theta) = \
  & = f_(Y_p,Y_(p - 1),...,Y_1) (y_p,...,y_(1);theta) times product^T_(t = p + 1) f_(Y_t vert Y_(t - 1),...,Y_(t - p)) (y_t vert y_(t - 1),...,y_(t - p);theta)
$
А функция правдоподобия
$
  ln L(theta vert Vect(y)) & = - T/2ln (2 pi) - T/2ln (sigma^2) + 1/2ln | Vect(V)^(- 1)_(p) | - \
  & wide - frac(1, 2 sigma^2) (Vect(y)_p - Vect(mu)_p)'Vect(V)^(- 1)_p (Vect(y)_p - Vect(mu)_p) - \
  & wide - sum^T_(t = p + 1) frac((y_t - alpha_0 - sum^p_(i = 1) alpha_i y_(t - i))^2, 2 sigma^2)
$
Если проигнорировать первые $p$ наблюдений, что допускается на больших выборках, то функция правдоподобия примет вид
$
  & ln L(theta vert y_T, ..., y_(p + 1)) = \
  & wide = - frac(T - p, 2) ln (2 pi) - frac(T - p, 2) ln (sigma^2) - sum^T_(t = p + 1) frac((y_t - alpha_0 - sum^p_(i = 1) alpha_i y_(t - i))^2, 2 sigma^2)
$
Эта функция дает нам в итоге МНК-оценки коэффициентов.
Оценка $sigma^2$
равна
$ hat(sigma)^2 = frac(1, T - p) sum^T_(t = p + 1) e^2_t $

== MA(1)

Рассмотрим процесс
$ Y_t = mu + epsilon_t + theta epsilon_(t - 1) $
где $epsilon sim "WN"(0, sigma^2)$.
Оцениваемый вектор параметров $theta = (mu,theta,sigma^2)$.

Если бы мы точно знали $epsilon_(t - 1)$, то
$
  Y_t vert epsilon_(t - 1) & sim N(mu + theta epsilon_(t - 1), sigma^2) \ f_(Y_t vert epsilon_(t - 1)) & = frac(1, sqrt(2 pi sigma^2)) exp [ - frac((y_t - mu - theta epsilon_(t - 1))^2, 2 sigma^2) ]
$

Далее следим за руками.
Предположим, что мы *знаем*, что $epsilon_0 = 0$.
Тогда
$ Y_1 vert epsilon_0 sim N(mu, sigma^2) $
Более того, мы *достоверно*#footnote[С точностью до *неизвестного* параметра $mu$.] знаем $epsilon_1 = y_1 - mu$; это позволяет нам записать
$ f_(Y_2 vert Y_1,epsilon_0) = frac(1, sqrt(2 pi sigma^2)) exp [ - frac((y_2 - mu - theta epsilon_1)^2, 2 sigma^2) ] $

Вообще, знание $epsilon_0 = 0$ позволяет нам достоверно рассчитать
${epsilon_t}$. Так как $epsilon_t = y_t - mu - theta epsilon_(t - 1)$,
то итоговая функция правдоподобия имеет вид
$ ln L(theta) = - T/2ln (2 pi) - T/2ln (sigma^2) - sum^T_(t = 1) frac(epsilon^2_t, 2 sigma^2) $

Кстати, именно на данном этапе проявляется роль обратимости MA-процесса: если $abs(theta) lt.double 1$, то эффект весьма произвольного выбора $epsilon_0$
крайне быстро затухает.

Если при численном решении задачи получается значение $hat(theta)$ больше 1, то оно должно быть отброшено, а процесс подбора начат заново с $1 slash hat(theta)$ в качестве начального значения.

С другой стороны мы можем записать точное значение функции правдоподобия *без* волюнтаристических предположений.
Ковариационная матрица выборки из процесса MA(1)
$ Vect(Omega) = Ex (Vect(Y) - Vect(mu)) (Vect(Y) - Vect(mu))' $
имеет вид
$
  Vect(Omega) = sigma^2 mat(
    delim: "[", (1 + theta^2), theta, 0, dots.h.c, 0 ;
    theta, (1 + theta^2), theta, dots.h.c, 0 ;
    0, theta, (1 + theta^2), dots.h.c, 0 ;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    0, 0, 0, dots.h.c, (1 + theta^2)
  )
$
А функция правдоподобия
$
  ln L(theta) = - T/2ln (2 pi) - 1/2ln | Vect(Omega) | - 1/2(Vect(y) - Vect(mu))'Vect(Omega)^(- 1) (Vect(y) - Vect(mu))
$ <eq-ma-mle>

=== Треугольное разложение

Гамильтон #cite(<hamilton1994time>) рассматривает следующую модификацию
@eq-ma-mle.
Известно#footnote[См. #cite(<hamilton1994time>, supplement: [стр.~87]).], что любая симметрическая положительно определенная матрица $Vect(Omega)$ единственным образом раскладывается следующим образом
$ Vect(Omega) = Vect(A) Vect(D) Vect(A)' $
где $Vect(A)$ — нижняя унитреугольная матрица#footnote[Треугольная матрица с единицами на главной диагонали.], $Vect(D)$ — диагональная матрица.

Из этого разложения следует, что
$ bar(Vect(Omega)) = bar(Vect(A)) times bar(Vect(D)) times bar(Vect(A)') = bar(Vect(D)) $
Также можно ввести следующие переменные
$ overline(Vect(y)) = Vect(A)^(- 1) (Vect(y) - Vect(mu)) $

С учетом этого @eq-ma-mle примет вид
$ ln L(theta) = - T/2ln (2 pi) - 1/2ln | Vect(D) | - 1/2overline(Vect(y))' Vect(D)^(- 1) overline(Vect(y)) $

Если обратить внимание на то, что
$
  Vect(A) Vect(D) Vect(A)'= underbrace(Vect(A) Vect(D)^(1 slash 2), Vect(L)) underbrace(Vect(D)^(1 slash 2) Vect(A)', Vect(L)')
$
где $Vect(L)$ — матрица разложения Холецкого, $Vect(D)^(1 slash 2)$ —
матрица, состоящая из корней элементов матрицы $Vect(D)$, и ввести переменную
$ tilde(Vect(y)) = Vect(L)^(- 1) (Vect(y) - Vect(mu)) $
то @eq-ma-mle примет вид
$ ln L(theta) = - T/2ln (2 pi) - ln | Vect(L) | - 1/2tilde(Vect(y))' tilde(Vect(y)) $

== ARMA(p, q)

Как ни странно, но ARMA(p, q) по построению функции правдоподобия ближе к процессам скользящего среднего.
Если взять начальные значения $Vect(y)_0 = (y_0,y_(- 1),...,y_(- p + 1))'$ и $Vect(epsilon)_0 = (epsilon_0,epsilon_(- 1),...,epsilon_(- q + 1))'$, то последовательность ${epsilon_t}$ может быть "достоверно" рассчитана
$ epsilon_t = y_t - alpha_0 - sum^p_(i = 1) alpha_i y_(t - i) - sum^q_(j = 1) theta_j epsilon_(t - j) $
что дает нам
$
  ln L(theta) & = ln f_(Y_T,Y_(T - 1),...,Y_1 vert Y_0,epsilon_0) (y_t,...,y_1 vert Vect(y)_0,Vect(epsilon)_(0);theta) = \
              & = - T/2ln (2 pi) - T/2ln (sigma^2) - sum^T_(t = 1) frac(epsilon^2_t, 2 sigma^2)
$

В качестве начальных значений можно взять:
- математические ожидания, то есть
  $
          y_t & = frac(alpha_0, 1 - sum^p_(i = 1) alpha_i), & t & = 0,- 1,...,- p + 1 \
    epsilon_t & = 0                                         & t & = 0,- 1,...,- q + 1
  $
- Бокс и Дженкинс #cite(<box2016time>) предложили брать в качестве $Vect(y)_0$ первые $p$ значений $y_t$, а в качестве $Vect(epsilon)_0$ — $epsilon_t = 0, t = p,...,p - q + 1$.
  Понятно, что расчет последовательности ${epsilon_t}$ будет проходить на наблюдениях $p + 1,...,T$; соответствующим образом меняется функция правдоподобия
  $
    ln L(theta) & = ln f_(Y_T,Y_(T - 1),...,Y_(p + 1) vert Y_0,epsilon_0) (y_t,...,y_(p + 1) vert Vect(y)_0,Vect(epsilon)_(0);theta) = \
    & = - frac(T - p, 2) ln (2 pi) - frac(T - p, 2) ln (sigma^2) - sum^T_(t = p + 1) frac(epsilon^2_t, 2 sigma^2)
  $

= GMM

Обобщенный метод моментов является, как понятно из названия, обобщением метода моментов.
Опишем этот метод согласно обзорной работе Хансена #cite(<hansen2007generalized>).

Рассмотрим наши данные, которые являются конечной реализацией процесса ${y_t}^T_(t = 1)$.
Сама модель определяется следующим простым
выражением
$ Ex f(y_(T); beta_0) = 0 $
где $beta_0 in PP subset RR^k$ — истинные значения параметров модели#footnote[Обычно к GMM прибегают в случае полупараметрических моделей, когда функцию правдоподобия записать "затруднительно".], а $f(...)$ — вектор-функция, имеющая $r$ компонент.

Для того, чтобы решение задачи было единственным, необходимо выполнение условия
$ Ex f(y_(t); beta) = 0 <=> beta = beta_0 $
$f(...)$ называют моментной функцией.
Реальное оценивание производится на выборочном аналоге матожидания
$ g_T (beta) = 1/T sum^T_(t = 1) f(y_(t);beta) $

В литературе показано, что
$ sqrt(T) g_T (beta_0) limits(->)^d N(0,V) $
где
$ V = sum^(infinity)_(j = - infinity) Ex [f(y_(t);beta_0) f(y_(t + j);beta_0)'] $
То есть $V$ — спектральная плотность на нулевой частоте.
Если
- $r < k$, то модель недоопределена; оценки в данном случае получить нельзя.
- $r = k$, то модель определена точно; у модели *может* быть решение.
  Это традиционный метод моментов.
- $r > k$, то модель переопределена; в общем случае решения также нет.
  Но возможно найти не точное, а максимально к нему приближенное решение.
  Это и есть обобщенный метод моментов.

== Оценивание модели

Одним из подходов к оцениванию является минимизация квадратичной формы
$ J(beta) = g_T (beta)'Vect(W) g_T (beta) $
то есть оценка равна
$ b_(italic("GMM")) = arg min_(beta in PP) J(beta) $
Данный подход обосновывается тем, что $J(beta)$ является выборочным
аналогом
$ Ex f(y_(t); beta)'Vect(W) Ex f(y_(t); beta) $

Полученная оценка является состоятельной; в работе #cite(<hansen1982large>)
было показано, что если процесс порождения данных стационарен и эргодичен,
то выполняется специальная форма закона больших чисел, по которой
$ sup_(beta in PP) | g_T (beta) - Ex f(y_(t); beta) | = 0 $
причем $PP$ компактно#footnote[Компактное пространство — топологическое пространство, в любом покрытии
  которого открытыми множествами найдётся конечное подпокрытие.].
Свойство компактности важно для "хорошего" поведения функции "вдали" от $beta_0$.

В качестве $Vect(W)$ часто предпочитают брать $Vect(V)^(- 1)$, так как в этом случае
$ T g_T (beta)'Vect(V)^(- 1) g_T (beta) limits(->)^d chi^2 (r) $

Оценкой ковариационной матрицы $b_(italic("GMM"))$ является #cite(<matyas1999generalized>, supplement: [гл.
  2])
$
  hat(Sigma)_(b_(italic("GMM"))) & = 1/T[F'_T Vect(V)^(- 1)_T F_T]^(- 1) \ F_T & = 1/T sum^T_(t = 1) frac(partial f(y_(t);beta), partial beta') bar.v_(beta = b_(italic("GMM")))
$
а $Vect(V)_T$ — оценка $Vect(V)$.

На практике для оценивания модели используется следующая последовательная процедура:

+ Задайте начальное значение $Vect(W)_0$. Удобно для этой цели взять единичную матрицу.
+ Минимизируйте $J(beta; Vect(W)_0)$
  $ b_1 = arg min_(beta in PP) g_T (beta)'Vect(W)_0 g_T (beta) $
+ Найдите новое значение $Vect(W)_1$
  $ Vect(W)_1 = { 1/T sum^T_(t = 1) g_T (b_1) g_T (b_1)' }^(- 1) $
+ Повторяйте шаги 2–3 до достижения сходимости.

== Примеры

=== Линейная регрессия

Рассмотрим пример линейной регрессии
$ y_t = beta_0 + beta_1 x_(1,t) + beta_2 x_(2,t) + ...+ beta_k x_(k,t) + epsilon_t $
где выполняются все предпосылки МНК и условия теоремы Гаусса-Маркова.

Из курса эконометрики известно, что для МНК оценок модели#footnote[В том числе *без* свободного члена.]
выполняется
$
  sum^T_(t = 1) x_(j,t) e_t = sum^T_(t = 1) x_(j,t) (y_t - hat(beta)_0 - sum^k_(i = 1) beta_i x_(i,t)) = 0, j = 1,..., k
$
Очевидно, что в качестве моментных функций можно взять
$ x_(j,t) (y_t - hat(beta)_0 - sum^k_(i = 1) beta_i x_(i,t)) $

=== Инструментальные переменные

Из теории метода инструментальных переменных известно, что инструменты должны удовлетворять двум свойствам:

+ релевантности, то есть $Cov(z_t, x_t) != 0$,
+ экзогенности, то есть $Cov(z_t, epsilon_t) = 0$.

Выборочный аналог второго свойства принимает вид
$ sum^T_(t = 1) z_t (y_t - beta_0 - beta_1 x_t) = 0 $
То есть в качестве моментных функций можно взять
$ z_t (y_t - beta_0 - beta_1 x_t) $

=== ARMA(2, 1)

Взято из блога #link("https://allenfrostline.com/blog/gmm-estimation/")[_Янга (Аллена) Вонга_].
Рассмотрим модель ARMA(2,1)
$ y_t = alpha_1 y_(t - 1) + alpha_2 y_(t - 2) + underbrace(epsilon_t + theta_1 epsilon_t, nu_t) $
где $epsilon_t sim "WN"(0, 1)$.

Данная модель дает нам следующие кандидаты в моментные функции:
+ $Ex (nu_t) = 0$;
+ $Var(nu_t) = (1 + theta^2_1)$;
+ $Cov(nu_t, nu_(t - 1)) = theta_1$;
+ $Cov(nu_t, nu_(t - j)) = 0, j >= 2$.

Если ограничиться этими условиями, то $f(y_(t);beta)$ можно записать в виде
$
  f(y_(t); beta) = mat(
    delim: "[", y_t - alpha_1 y_(t - 1) - alpha_2 y_(t - 2);
    (y_t - alpha_1 y_(t - 1) - alpha_2 y_(t - 2))^2 - (1 + theta^2_1);
    (y_t - alpha_1 y_(t - 1) - alpha_2 y_(t - 2)) (y_(t - 1) - alpha_1 y_(t - 2) - alpha_2 y_(t - 3)) - theta_(1);
    (y_t - alpha_1 y_(t - 1) - alpha_2 y_(t - 2)) (y_(t - 2) - alpha_1 y_(t - 3) - alpha_2 y_(t - 4))
  )
$
Откуда
$ Ex f(y_t, hat(beta)) = 0 $

Для расчета ковариационной матрицы и для целей процедуры оценивания рассчитаем частные производные $f(y_(t);beta)$:
$ frac(partial f, partial theta') = mat(delim: "[", cal(E)_1, cal(E)_2, cal(E)_3) $
$
  cal(E)_1 = mat(
    delim: "[", - y_(t - 1);
    - 2 y_(t - 1) (y_t - alpha_1 y_(t - 1) - alpha_2 x_(t - 2));
    - y_(t - 1) (y_(t - 1) - alpha_1 y_(t - 2) - alpha_2 y_(t - 3)) - y_(t - 2) (y_t - alpha_1 y_(t - 1) - alpha_2 y_(t - 2));
    - y_(t - 1) (y_(t - 2) - alpha_1 y_(t - 3) - alpha_2 y_(t - 4)) - y_(t - 3) (y_t - alpha_1 y_(t - 1) - alpha_2 y_(t - 2))
  )
$
$
  cal(E)_2 = mat(
    delim: "[", - y_(t - 2);
    - 2 y_(t - 2) (y_t - alpha_1 y_(t - 1) - alpha_2 y_(t - 2));
    - y_(t - 2) (y_(t - 1) - alpha_1 y_(t - 2) - alpha_2 y_(t - 3)) - y_(t - 3) (y_t - alpha_1 y_(t - 1) - alpha_2 y_(t - 2));
    - y_(t - 2) (y_(t - 2) - alpha_1 y_(t - 3) - alpha_2 y_(t - 4)) - y_(t - 4) (y_t - alpha_1 y_(t - 1) - alpha_2 y_(t - 2))
  )
$
$
  cal(E)_3 = mat(
    delim: "[", 0 ;
    - 2 theta_(1);
    - 1;
    0
  )
$

#bibliography("bibliography.bib")

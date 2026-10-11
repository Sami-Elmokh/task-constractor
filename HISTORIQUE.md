# Où ChatGPT se trompe

ChatGPT n'a reçu que `problem.md`. Son code a ensuite été lancé contre `solution_test.py`.

## 1. `driven_oscillator_resonance` : ChatGPT échoue à 6/12

**Le problème.** Un oscillateur de fréquence propre $\omega$ quelconque, forcé par une force périodique donnée par $N$ échantillons :

$$y'' + \omega^2 y = \mu(t), \qquad y(0) = y'(0) = 0.$$

$\omega$ peut être exactement égal à un harmonique $k$ de la force, ou à $10^{-12}$ près. La précision demandée est $10^{-9}$, jusqu'à $t = 200$.

| Erreur de ChatGPT | Tests échoués |
|---|---|
| **Formule fausse à la résonance** : $C_k = \frac{t \sin kt}{2k} + \frac{\cos kt - 1}{2k^2}$. Le second terme ne résout pas l'équation ($C'' + k^2C - \cos kt = -\tfrac12$ au lieu de 0). | 3 (résonance exacte) |
| **Bascule sur la formule de résonance** dès que $\lvert\omega - k\rvert < 10^{-7}$. Près de la résonance, ce n'est pas assez précis : même avec la bonne formule, ces tests échouent. | 3 ($\omega$ à $10^{-12}$ et à $5\times10^{-9}$ d'un harmonique) |

## 2. `oscillator_period_energy` : ChatGPT échoue à 5/14

**Le problème.** Une particule dans un puits de potentiel $f(x)$, avec $x'' = -f'(x)$. Il faut calculer la période en fonction de l'énergie, des petites oscillations jusqu'au seuil d'échappement :

$$T(h) = 2\int_a^b \frac{du}{\sqrt{2\,(h - f(u))}}.$$

La précision demandée est $10^{-8}$, y compris pour $h$ à $10^{-6}$ du sommet d'une barrière et pour $h - f(x_0)$ jusqu'à $10^{-8}$.

| Erreur de ChatGPT | Tests échoués |
|---|---|
| **Saute par-dessus une barrière étroite** : son pas de recherche grandit, donc il franchit la barrière du pendule, large de seulement $2{,}8\times10^{-3}$. Il obtient des points de rebroussement à $\pm 2{,}3\times10^{17}$ au lieu de $\pm 3{,}14$, et une période de $10^{18}$ au lieu de 34,56. Sur le potentiel incliné, le code plante. | 3 |
| **Annulation numérique** : il calcule $h - f(u)$ directement. Pour les toutes petites oscillations, l'erreur est de $5\times10^{-5}$ alors que la précision demandée est $10^{-8}$. | 2 |

## 3. `rectangular_pair_normal_form` : ChatGPT échoue à 8/14

**Le problème.** Un couple de matrices entières $A$ ($n\times m$) et $B$ ($m\times n$), donné à équivalence simultanée près ($A \mapsto QAP^{-1}$, $B \mapsto PBQ^{-1}$). Il faut calculer sa forme normale : les chaînes nilpotentes, par longueur et par côté de départ, et les blocs de Jordan $(\lambda, r)$ de $BA$ pour $\lambda \neq 0$, puis le couple canonique. Les coefficients atteignent environ $10^8$, donc les rangs doivent être exacts.

| Erreur de ChatGPT | Tests échoués |
|---|---|
| **Formule des chaînes fausse** : il lit le nombre de chaînes de chaque côté dans les différences premières des chutes de rang de ce côté, sans tenir compte de l'alternance entre les deux côtés. | 2 |
| **Confusion entre $m$ et $n$** : $BA$ est de taille $m\times m$, mais il utilise $n$. Le code plante dès que $m \neq n$. | 6 (blocs inversibles et forme canonique) |

Toutes ces erreurs portent sur des régimes **écrits dans `problem.md`**. Les échecs viennent donc de vraies erreurs, pas d'une ambiguïté de l'énoncé.

# Les erreurs de ChatGPT sur `rademacher_sum_moments`

ChatGPT n'a reçu que `problem.md`. Son code échoue à **6 des 16 tests**, pour deux erreurs indépendantes.

| Fonction | Tests réussis | Cause des échecs |
|---|---|---|
| `rademacher_law` | 5/5 | — |
| `even_moment` | 0/4 | Erreur 1 : récurrence moments–cumulants fausse |
| `tail_probability` | 3/3 | — |
| `khintchine_ratio` | 2/4 | Erreur 1 ($p = 2$ et $p = 4$), puis erreur 2 ($p = 40$) |

## Erreur 1 : la récurrence moments–cumulants ne garde que les ordres impairs

**Rappel.** Si $K(t) = \log E[e^{tS}] = \sum_j \kappa_j t^j / j!$ (les cumulants $\kappa_j$), les moments $\mu_r = E[S^r]$ vérifient

$$\mu_r = \sum_{j=1}^{r} \binom{r-1}{j-1} \kappa_j\, \mu_{r-j}.$$

Ici $E[e^{tS}] = \prod_i \cosh(c_i t)$ est une fonction **paire** : $S$ est symétrique, donc **tous les cumulants d'ordre impair sont nuls**. Seuls les $\kappa_j$ avec $j$ pair contribuent.

**Ce que fait ChatGPT.** Ses cumulants pairs sont justes : $\kappa_{2k} = \frac{2^{2k}(2^{2k}-1)B_{2k}}{2k}\sum_i c_i^{2k}$, ce qui donne bien $\kappa_2 = \sum c_i^2$ et $\kappa_4 = -2\sum c_i^4$. Mais dans la récurrence, il fait la somme sur les ordres **impairs** :

```python
for j in range(1, r + 1, 2):   # j = 1, 3, 5, ...  au lieu de 2, 4, 6, ...
```

Il ne garde donc que des $\kappa_j$ nuls, et $\mu_r = 0$ pour tout $r \ge 1$. **Tous ses moments pairs valent 0**, ce qui est absurde puisque $E[S^{2m}] \ge E[S^2]^m > 0$.

**Contre-exemple minimal : $c = (1, 2, 3)$, $m = 2$.** $S$ vaut $0$ avec probabilité $1/4$ et $\pm 2, \pm 4, \pm 6$ avec probabilité $1/8$ chacun, donc

$$E[S^4] = \frac{2(2^4 + 4^4 + 6^4)}{8} = 392 = 3\cdot 14^2 - 2\cdot 98.$$

ChatGPT renvoie $0$.

**Résultat sur les tests :**

| Test | ChatGPT | Attendu |
|---|---|---|
| `even_moment` normal, $c = (1,2,3)$, $m = 2$ | 0 | 392 |
| `even_moment` limite, $c = (7,-3,0,2)$, $m = 1$ | 0 | 62 |
| `even_moment` 1000 coefficients égaux, $m = 10$ | 0 | 635352202681178618159158926784135168000 |
| `even_moment` $10^5$ coefficients, $m = 2$ et $m = 20$ | 0 | (entiers exacts de 34 et 354 chiffres) |
| `khintchine_ratio` $p = 2$ | 0.0 | 1.0 |
| `khintchine_ratio` $p = 4$, 1000 coefficients égaux | 0.0 | 1.3158546124259147 |

**Remarque.** Son code calcule aussi les coefficients de $\tanh$ et une liste `log_coeff`, qu'il **n'utilise jamais**. Le commentaire annonce une méthode, puis le code en suit une autre.

## Erreur 2 : conversion en `float` d'un moment hors de portée

**Ce que fait ChatGPT.** Pour $p$ pair, il calcule le rapport ainsi :

```python
return float(moment) ** (1.0 / p) / variance ** 0.5
```

**Conséquence.** Pour $10^5$ coefficients et $p = 40$, $E[S^{40}] \approx 9 \times 10^{353}$, bien au-delà du plus grand flottant ($\approx 1{,}8 \times 10^{308}$) : `OverflowError: int too large to convert to float`. Le rapport cherché vaut pourtant seulement $3{,}87$.

**Vérification.** En corrigeant uniquement la boucle de l'erreur 1 (`range(2, r + 1, 2)`), son code passe 15 tests sur 16. Le seul échec restant est ce test $p = 40$. Les deux erreurs sont donc bien indépendantes.

**Ironie.** Pour $p$ non pair, il travaille en logarithmes justement « pour éviter le dépassement ». Il n'applique pas cette précaution à la branche $p$ pair, la seule où les nombres deviennent énormes.

## Ce que ChatGPT a bien fait

- **Loi exacte.** Il a calculé la loi par programmation dynamique avec des comptes entiers exacts. Les coefficients nuls, les signes opposés et les valeurs qui coïncident sont bien traités : 5 tests sur 5.
- **Queue stricte.** Il utilise bien $\lvert S\rvert > t$ et pas $\ge t$, même quand $t$ est exactement une valeur prise par $\lvert S\rvert$. Il convertit aussi $t$ en fraction exacte : 3 tests sur 3.
- **Normes impaires et non entières.** Il calcule $E\lvert S\rvert^p$ à partir de la loi, et non $E[S^p]$, pour $p = 3$, $7$ et $1{,}5$.

## Pourquoi ces échecs sont légitimes

- **Les exigences sont dans `problem.md`.** Les moments doivent être exacts, avec jusqu'à $10^5$ coefficients et $m \le 20$, et le rapport doit être calculé pour $p$ pair sur ces mêmes entrées. L'énoncé indique aussi que les $\varepsilon_i$ sont des signes symétriques.
- **L'erreur 1 est mathématique.** Une variable symétrique n'a que des cumulants pairs. Garder les impairs donne des moments nuls, ce qui est faux pour toute entrée non nulle, même la plus simple.
- **L'erreur 2 est une vraie erreur de conception, pas un manque de précision.** Le résultat n'est pas légèrement faux : le code plante, alors que le rapport cherché vaut $3{,}87$ et se calcule exactement en fraction avant une seule conversion en `float`.

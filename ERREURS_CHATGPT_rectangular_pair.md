# Les erreurs de ChatGPT sur `rectangular_pair_normal_form`

ChatGPT n'a reçu que `problem.md`. Son code échoue à **8 des 14 tests**, pour deux erreurs indépendantes.

| Fonction | Tests réussis | Cause des échecs |
|---|---|---|
| `alternating_ranks` | 3/3 | — (rangs exacts, aucune erreur) |
| `chain_counts` | 2/4 | Erreur 1 : formule des chaînes fausse |
| `invertible_blocks` | 1/4 | Erreur 2 : confusion entre $m$ et $n$ |
| `canonical_form` | 0/3 | Conséquence des erreurs 1 et 2 |

## Erreur 1 : la formule des chaînes ne tient pas compte de l'alternance

**Rappel.** Une chaîne de longueur $L$ est une suite de vecteurs $e_0 \to e_1 \to \dots \to e_{L-1} \to 0$ qui **alternent** entre $\mathbb{C}^m$ et $\mathbb{C}^n$ : on passe d'un côté à l'autre par $A$, puis on revient par $B$.

**Ce que fait ChatGPT.** Il traite chaque côté séparément. Pour le côté $m$, il prend les chutes de rang $d_k = \mathrm{rank} W^m_k - \mathrm{rank} W^m_{k+1}$ et pose « nombre de chaînes de longueur $L$ partant de $m$ $= d_{L-1} - d_L$ ». Il fait de même pour le côté $n$. C'est la formule du cas d'**une seule** matrice nilpotente, comme les blocs de Jordan de la question 10, mais elle ne s'applique pas ici.

**Contre-exemple minimal : une seule chaîne de longueur 3 partant de $m$.** On a $e_0 \in \mathbb{C}^m$, $e_1 = Ae_0 \in \mathbb{C}^n$, $e_2 = Be_1 \in \mathbb{C}^m$, puis $Ae_2 = 0$.

| $k$ | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| $\mathrm{rank} W^m_k$ ($I$, $A$, $BA$, $ABA$) | 2 | 1 | 1 | 0 |
| $\mathrm{rank} W^n_k$ ($I$, $B$, $AB$, $BAB$) | 1 | 1 | 0 | 0 |
| **total** | **3** | **2** | **1** | **0** |

- La suite d'**un seul côté** oscille (2, 1, 1, 0), parce que les vecteurs de la chaîne changent de côté à chaque pas. Les différences de ses chutes ne comptent donc rien de sensé. ChatGPT trouve `{(1,'m'): 1, (2,'n'): 1}`, c'est-à-dire deux chaînes qui n'existent pas, au lieu de `{(3,'m'): 1}`.
- La bonne méthode **additionne les deux côtés**. Le total (3, 2, 1, 0) vaut exactement $\max(L - k, 0)$, et sa différence seconde donne bien une chaîne de longueur 3. Le côté de départ se lit ensuite dans la **différence** des deux côtés, où les blocs inversibles s'annulent.

**Un second défaut dans le même code.** Il ne cherche que des longueurs $L \le \max(m, n)$, alors qu'une chaîne peut être plus longue. Dans l'exemple, $L = 3 > \max(2, 1)$, donc la vraie chaîne ne pouvait même pas être trouvée.

**Résultat sur les tests :**

| Test | ChatGPT | Attendu |
|---|---|---|
| normal | `{(1,'m'): 2, (2,'n'): 2, (3,'m'): 1}` | `{(3,'m'): 1, (2,'n'): 1}` |
| difficile | `{(1,'n'): 4, (2,'m'): 3, (3,'n'): 2, (4,'m'): 1, (5,'m'): 1, (5,'n'): 1}` | `{(1,'m'): 1, (1,'n'): 1, (2,'m'): 1, (3,'n'): 1, (4,'m'): 1, (5,'m'): 1, (5,'n'): 1}` |

**Remarque.** Son code construit ensuite un système d'équations « pour vérifier la cohérence », mais **ne le résout jamais** et renvoie le résultat faux. Le commentaire annonce une vérification qui n'a pas lieu.

## Erreur 2 : confusion entre $m$ et $n$ dans `invertible_blocks`

**Rappel.** $A$ est de taille $n \times m$ et $B$ de taille $m \times n$, donc $BA$ est de taille $m \times m$.

**Ce que fait ChatGPT.** Il pose `n = len(A)`, le nombre de lignes de $A$, donc $n$, puis l'utilise comme taille de $BA$ : la matrice identité, le polynôme caractéristique et les puissances de $BA - \lambda I$ sont tous dimensionnés avec $n$.

**Conséquence.** Dès que $m \neq n$, le code plante (`IndexError: list index out of range`). Dans le test normal, $m = 4$ et $n = 3$. Le seul test réussi est celui où $m = n = 4$.

**Vérification.** En remplaçant seulement $n$ par $m$ dans cette fonction, son code donne les bons blocs sur les tests normal et difficile. Sa méthode (polynôme caractéristique exact, puis rangs de $(BA - \lambda I)^k$) était donc juste : c'est une erreur de programmation sur des dimensions clairement indiquées dans `problem.md`.

## Ce que ChatGPT a bien fait

- **Rangs exacts.** Il a calculé les rangs en arithmétique entière exacte, pas en flottants, donc il a évité le piège des coefficients d'environ $10^8$, où `numpy.linalg.matrix_rank` se trompe.
- **Forme canonique.** Son assemblage de $(A_0, B_0)$ respecte l'ordre et la disposition demandés. Ce sont les erreurs 1 et 2, en amont, qui le font échouer.

## Pourquoi ces échecs sont légitimes

- **Les définitions sont dans `problem.md`.** Les chaînes y sont définies explicitement comme alternant entre $\mathbb{C}^m$ et $\mathbb{C}^n$, et les formes $n\times m$ et $m\times n$ sont écrites.
- **Ce sont deux erreurs réelles.** L'erreur 1 est mathématique : elle applique la formule d'un nilpotent seul à un couple d'applications alternées. L'erreur 2 est une erreur de code. Aucune ne vient d'une ambiguïté de l'énoncé.
- **L'erreur 1 correspond à la difficulté centrale du sujet.** C'est exactement la combinatoire des questions 14c et 17 de l'épreuve X-ENS 2025.

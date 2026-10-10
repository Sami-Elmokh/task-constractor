# Historique des modifications, étape par étape

Ce document reprend **chaque commit** du dépôt dans l'ordre chronologique. Pour chaque étape, il explique :

- **ce qui a changé** ;
- **pourquoi** ;
- **comment c'est vérifié** ;
- **la commande pour voir le changement exact**, à lancer depuis le dossier du dépôt.

Pour suivre en parallèle, `git log --oneline` liste les commits et `git show <hash>` affiche le contenu exact de l'un d'eux.

| # | Commit | Date | Résumé |
|---|---|---|---|
| 1 | `2cb0c4f` | 05/10 | Première version : deux tâches au format SciCode (SME01, SME02) |
| 2 | `b6b3c00` | 05/10 | Passage au format du guide SME (6 fichiers par tâche) |
| 3 | `4d5ef96` | 05/10 | README mis à jour |
| 4 | `20d382a` | 05/10 | Alignement exact sur le template |
| 5 | `ec75d3a` | 05/10 | Ta réécriture de `verification.md` (tâche 2) |
| 6 | `650d27e` | 06/10 | Tests de régimes rétablis et vocabulaire « deliberate error » |
| 7 | `546fba0` | 06/10 | `problem.md` réduit aux exigences, sans indices |
| 8 | `9d6a035` | 10/10 | Durcissement de la tâche 1 après le retour du relecteur |
| 9 | `05a2aa9` | 10/10 | Contrôle par modèle de pointe de la tâche 2 enregistré |

---

## Étape 1 : première version au format SciCode (`2cb0c4f`)

**Ce qui a changé.** Deux tâches ont été créées, au format du dépôt officiel SciCode (`problem.json`, `test_data.h5`, dossier `gold/`, scripts `build_task.py`, `run_tests.py` et `validate.py`, plus un `TASK.md`) :

- **SME01** : oscillateur forcé $y'' + y = \mu(t)$, réponse périodique ou résonance ;
- **SME02** : période d'un oscillateur conservatif $x'' = -f'(x)$.

**Pourquoi.** À ce moment-là, le guide officiel n'était pas encore disponible. On a donc repris le format de SciCode lui-même, découvert en travaillant sur la tâche 9.1 (Weighted Jacobi).

**Ce qu'on a appris.** En vérifiant la tâche 9.1, on a découvert que **ses réponses attendues sont fausses** pour $\omega \ne 1$ : elles viennent d'une itération qui oublie le terme $(1-\omega)\,x$. Ça a montré **pourquoi une solution de référence doit être vérifiée par des méthodes indépendantes**. Ce principe a guidé toute la suite.

**Bugs trouvés pendant cette étape, grâce aux vérifications indépendantes :**

- `np.fft.fftfreq(N, 1/N)` renvoie parfois `1.0000000000000002` au lieu de `1`, par exemple pour $N = 49$. Le mode résonant n'était alors pas détecté, et le code divisait par environ $10^{-16}$.
- Une recherche des points de rebroussement à pas fixe **sautait par-dessus une barrière étroite** et donnait une période de 241 au lieu de 34,56.

**Voir le changement :** `git show 2cb0c4f --stat`

---

## Étape 2 : passage au format du guide SME (`b6b3c00`)

**Ce qui a changé.** Le guide (`scicode_sme_task_guide.md`) et le template sont arrivés. Les deux tâches ont été **entièrement réécrites** au format demandé, avec 6 fichiers par tâche :

| Fichier | Rôle |
|---|---|
| `problem.md` | L'énoncé complet, sans solution ni réponses |
| `solution.py` | Le code complet |
| `main.py` | Lance la solution sur un exemple |
| `solution_test.py` | Les tests `pytest`, avec les réponses attendues **écrites en dur** |
| `requirements.txt` | Les bibliothèques utilisées |
| `verification.md` | Comment la science et les réponses ont été vérifiées |

Les anciens fichiers (`problem.json`, `.h5`, `gold/`, scripts) ont été supprimés, comme tu l'avais demandé.

**Pourquoi chaque choix :**

- **Originalité.** Le guide interdit de « copier ou réécrire légèrement un exercice en ligne ». Les exemples qui venaient des sujets de concours, comme $|\sin t|$ ou $-\cos x + x^2/20$, ont donc été **remplacés par nos propres exemples** : $e^{\cos(t-0.4)}$, un pendule incliné, un puits quartique et un réseau dans un piège.
- **La dernière fonction doit utiliser les précédentes**, comme le guide l'exige. Pour la tâche 2, on a ajouté `period_energy_curve`, qui appelle `barrier_energy` et `oscillation_period`.
- **Les réponses attendues ne sont jamais calculées avec `solution.py`**, autre exigence du guide. Elles viennent de scripts indépendants : formules exactes, séries de Bessel, intégrales elliptiques, intégration directe de l'EDO.

**Bug trouvé pendant cette étape.** Le premier « réseau piégé » ($0{,}05x^4 - \cos 2x$) **n'avait aucune barrière intérieure**, parce que le terme $x^4$ les effaçait. Le test qui prétendait que « l'orbite traverse trois puits » était donc faux. Il a été remplacé par $0{,}01x^4 - \cos 2x$, dont on a vérifié numériquement les barrières, de hauteur 1,064 en $x = \pm 1{,}613$.

**Comment c'est vérifié.** Pour chaque tâche : `python main.py`, `python -m pytest -q`, et le contrôle d'erreurs volontaires décrit à l'étape 6.

**Voir le changement :** `git show b6b3c00 --stat`

---

## Étape 3 : README (`4d5ef96`)

**Ce qui a changé.** Le `README.md` décrit maintenant les deux nouvelles tâches et les 6 fichiers du format SME.

**Pourquoi.** Au commit précédent, la mise à jour du README n'avait pas été enregistrée, parce que le fichier avait changé sur le disque pendant l'opération. Le README poussé décrivait donc encore l'ancien format, d'où cette correction.

**Voir le changement :** `git show 4d5ef96`

---

## Étape 4 : alignement exact sur le template (`20d382a`)

**Ce qui a changé.** Chaque fichier a été aligné sur la **structure** du template, en gardant des noms de fonctions parlants comme tu l'as demandé :

- `solution.py` : même docstring que le template, plus de fonction d'aide ni de constante en dehors des fonctions demandées ;
- `solution_test.py` :
  - même docstring et constantes `RTOL` / `ATOL` en tête ;
  - chaque test suit le schéma `test_input` → `expected` → `actual` → `assert` ;
  - les tests de la fonction finale s'appellent `test_complete_solution_*` ;
- `main.py` : `example_input` → `result` → `print(result)` ;
- `verification.md` : les 4 puces du template.

**Pourquoi.** Tu as demandé que le résultat final ait **exactement** le même format que le template, en nombre de fichiers et en structure.

**Voir le changement :** `git show 20d382a`

---

## Étape 5 : ta réécriture de `verification.md` (`ec75d3a`)

**Ce qui a changé.** Tu as restructuré `oscillator_period_energy/verification.md` avec des sous-titres en gras et des listes, ce qui le rend plus lisible.

**Ce qu'il faut retenir.** À l'étape suivante, mon push a été **refusé** par GitHub, parce que ton commit n'était pas encore dans ma copie locale. Je n'ai rien forcé : j'ai placé mon commit **par-dessus le tien** (`git rebase`), j'ai gardé ta version du fichier et j'y ai seulement reporté mes changements de contenu.

**Voir le changement :** `git show ec75d3a`

---

## Étape 6 : tests de régimes rétablis et vocabulaire (`650d27e`)

**Ce qui a changé :**

1. **Deux tests rétablis** dans la tâche 2 :
   - `test_turning_points_asymmetric_well_case` : un puits incliné, donc asymétrique ($a \ne -b$) ;
   - `test_oscillation_period_multiwell_case` : une orbite qui passe au-dessus de deux barrières intérieures.
2. Dans `verification.md`, « injected errors » est devenu **« Deliberate-error check »**.

**Pourquoi :**

1. À l'étape 4, j'avais réduit à exactement 3 tests par fonction. Mais le guide dit : *« Add more when the function has other scientific regimes or branches »*. Ces deux tests couvrent des régimes que les autres n'atteignent pas.
2. Le guide parle de *« catch a deliberate error such as a wrong sign or skipped step »*. Le mot « injection » prêtait à confusion.

**Ce qu'est le « contrôle d'erreurs volontaires ».** On copie la tâche, on introduit **une seule** erreur réaliste dans la copie de `solution.py` (un signe oublié, une demi-période, une formule approchée…), puis on lance les tests. Si au moins un test échoue, les tests attrapent bien cette erreur. C'est fait pour chaque erreur, l'une après l'autre.

**Voir le changement :** `git show 650d27e`

---

## Étape 7 : `problem.md` réduit aux exigences (`546fba0`)

**Ce qui a changé.** Les deux `problem.md` ont été raccourcis : de 851 à 571 mots pour la tâche 2, de 620 à 505 mots pour la tâche 1. On a retiré tout ce qui expliquait **comment** contourner les difficultés.

**Pourquoi.** Un énoncé trop détaillé **rend la tâche facile** : il donne au modèle la liste des pièges. Le guide demande seulement de dire **quoi** calculer.

| On garde (exigé par le guide) | On retire (des indices) |
|---|---|
| Les définitions (intervalle maximal, interpolant…) | « la barrière peut faire $10^{-3}$ de large » |
| Les formules qui définissent le calcul | « singularité intégrable en $1/\sqrt{\cdot}$ » |
| Les hypothèses | « attention à l'annulation numérique » |
| Les régimes couverts et la précision exigée | Les formules explicites de $a_1$ et $b_1$, qui donnaient la fonction 1 |

**La règle à retenir.** Chaque régime testé doit être **mentionné** dans `problem.md`, avec la précision exigée, sinon le test serait injuste. Mais il ne faut pas dire **comment** le traiter.

**Comment c'est vérifié.** Les tests n'ont pas changé, et on a contrôlé que chacun repose sur une exigence écrite dans l'énoncé.

**Voir le changement :** `git show 546fba0 -- '*/problem.md'`

---

## Étape 8 : durcissement de la tâche 1 (`9d6a035`)

**Le point de départ.** Le relecteur a répondu : *la tâche est correcte, mais un modèle de pointe l'a résolue et passe les 9 tests*. Il a suggéré trois pistes :

- une fréquence propre variable ;
- des tests à la résonance exacte et tout près de la résonance ;
- une exigence de précision près de la résonance.

**Ce qui a changé.**

1. **Fréquence propre variable :** l'équation devient

   $$y'' + \omega^2 y = \mu(t), \qquad y(0) = y'(0) = 0,$$

   avec $\omega > 0$ quelconque. La résonance peut maintenant se produire avec **n'importe quel harmonique** $k$ ($\omega = k$), plus seulement avec $k = 1$.
2. **De nouvelles fonctions :**

   | Fonction | Rôle |
   |---|---|
   | `forcing_coefficients(mu)` | Tous les coefficients $a_k$ et $b_k$ du forçage |
   | `mode_response(k, omega, t)` | La réponse depuis le repos à $\cos kt$ et à $\sin kt$ |
   | `driven_response(mu, omega, t)` | La somme sur tous les harmoniques |
3. **Des tests près de la résonance** : $\omega = k$ exactement, $\omega$ à $10^{-12}$ d'un harmonique (au-dessus et en dessous), et $\omega$ à $5\times10^{-9}$. La précision exigée est $10^{-9}$, jusqu'à $t = 200$.

**Pourquoi c'est difficile.** Il y a deux pièges.

- **La formule du cours est instable.** Pour $\omega \ne k$,

  $$C_k(t) = \frac{\cos kt - \cos \omega t}{\omega^2 - k^2}.$$

  Quand $\omega$ est très proche de $k$, le numérateur et le dénominateur sont tous les deux minuscules. Le numérateur est une différence de deux nombres presque égaux, donc on perd presque tous les chiffres significatifs (c'est l'**annulation numérique**). On a mesuré une erreur relative de $3\times10^{-4}$ pour $|\omega - k| = 10^{-12}$.
- **Basculer sur la formule résonante sous un seuil est faux aussi.** Si $|\omega - k| < \varepsilon$, utiliser $t \sin kt / (2k)$ donne une erreur de l'ordre de $|\omega - k|\,t$. Les deux tests (à $10^{-12}$ et à $5\times10^{-9}$) sont choisis pour qu'**aucun seuil** ne passe les deux.

**La solution stable.** On réécrit la formule pour **ne jamais diviser par $\omega - k$** :

$$\cos kt - \cos\omega t = 2\sin\tfrac{(\omega+k)t}{2}\,\sin\tfrac{(\omega-k)t}{2},$$

et on calcule $\sin(\delta t/2)/\delta$ avec une fonction sinc, qui vaut exactement $t/2$ quand $\delta = 0$. C'est exact à la résonance, et précis à $10^{-13}$ tout près d'elle. Cette astuce **n'est pas** donnée dans `problem.md`.

**Comment c'est vérifié :**

- **Réponses attendues :** calculées avec **50 chiffres significatifs** (bibliothèque `mpmath`), ce qui rend la formule du cours exacte, puis vérifiées par une intégration directe de l'EDO, avec un accord à $10^{-11}$.
- **Contrôle d'erreurs volontaires :** 11 erreurs réalistes, toutes détectées. Parmi elles, la formule du cours, le basculement pour chaque seuil de $10^{-6}$ à $10^{-12}$, et une résonance gérée seulement pour $k = 1$.
- **Contrôle par modèle de pointe :** ChatGPT, qui n'avait reçu que `problem.md`, **échoue à 6 des 12 tests** :
  - sa formule à la résonance est fausse : on obtient $C'' + k^2 C - \cos kt = -1/2$ au lieu de 0 ;
  - il bascule sur la formule résonante quand $|\omega - k| < 10^{-7}$. Même après correction de sa formule, il échoue encore à 3 tests.

**Voir le changement :** `git show 9d6a035 -- driven_oscillator_resonance/`

---

## Étape 9 : contrôle par modèle de pointe de la tâche 2 (`05a2aa9`)

**Ce qui a changé.** Le résultat du contrôle est ajouté à `oscillator_period_energy/verification.md`, dans ta mise en forme.

**Le résultat.** ChatGPT, qui n'avait reçu que `problem.md`, **échoue à 5 des 14 tests** :

| Erreur de ChatGPT | Conséquence |
|---|---|
| Son pas de recherche grandit géométriquement et **saute par-dessus la barrière étroite** du pendule | Points de rebroussement à $\pm 2{,}3\times10^{17}$ au lieu de $\pm 3{,}14$ ; période de $10^{18}$ au lieu de 34,56 ; plantage sur le potentiel incliné |
| Il calcule $h - f(u)$ directement | Pour les toutes petites oscillations, une erreur de $5\times10^{-5}$ alors que la précision exigée est $10^{-8}$ |

Ces deux régimes sont **explicitement mentionnés** dans `problem.md` : $h$ à $10^{-6}$ d'un sommet de barrière, et $h - f(x_0)$ jusqu'à $10^{-8}$. Les échecs viennent donc de vraies erreurs, pas d'une ambiguïté.

**Voir le changement :** `git show 05a2aa9`

---

## Les leçons générales

1. **Vérifier la solution de référence avec des méthodes indépendantes.** Elle peut être fausse, comme la tâche 9.1 de SciCode, ou mon bug `fftfreq`.
2. **Vérifier aussi les outils de vérification.** Mon premier script de contrôle était faux (il manquait `np`), et un de mes potentiels de test n'avait pas les barrières annoncées.
3. **Ne jamais affirmer sans avoir mesuré.** J'avais écrit qu'un `quad` simple passait tous les tests. C'était faux, et j'ai corrigé après l'avoir mesuré.
4. **Un énoncé dit *quoi*, pas *comment*.** Les régimes et la précision doivent être écrits ; les astuces, non.
5. **Une tâche est valide quand un modèle de pointe échoue pour une vraie raison scientifique ou numérique**, pas à cause d'une ambiguïté.
6. **Ne jamais forcer un push.** Quand GitHub refuse, on récupère les changements de l'autre personne et on les fusionne.

## État final

| | `driven_oscillator_resonance` | `oscillator_period_energy` |
|---|---|---|
| Tests de notre solution | 12/12 ✅ | 14/14 ✅ |
| Erreurs volontaires détectées | 11/11 | 9/9 |
| ChatGPT (avec seulement `problem.md`) | échoue à 6/12 → **valide** | échoue à 5/14 → **valide** |

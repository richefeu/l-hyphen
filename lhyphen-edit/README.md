lhedit
======

Un petit éditeur en terminal pour les fichiers d'entrée de l-hyphen (`input.txt`,
fichiers `conf*`) : coloration syntaxique, documentation en ligne des mots-clés et
snippets, assez léger pour être utilisé par ssh sur un cluster.

Aucune dépendance : termios et séquences ANSI uniquement, un compilateur C++17
suffit.

Compilation
-----------

~~~bash
cd lhyphen-edit
make                                   # compile lhedit et le copie à la racine (à côté de run et see2)
../lhedit ../examples/Poutre/input.txt
~~~

Touches
-------

Façon nano ; la barre du bas les rappelle, et `^G` ouvre la liste complète.

| Touche | Action |
|---|---|
| `^O` (ou `^S`) | enregistrer ; demande un nom s'il n'y en a pas |
| `^X` | quitter (propose d'enregistrer si le fichier est modifié) |
| `^W` / `^N` | chercher / occurrence suivante, sans tenir compte de la casse |
| `^L` | aller à une ligne |
| `^D` | afficher/masquer la documentation de la sélection, ou du mot sous le curseur |
| `^P` | insérer un snippet (taper pour filtrer la liste) |
| `Shift` + flèches | sélectionner ; aussi avec `Home`, `End`, `PgUp`, `PgDn` |
| `^K` / `^U` | couper la sélection, ou la ligne entière s'il n'y en a pas / coller |
| `^Z` / `^Y` | annuler / refaire |
| `^A` / `^E` | début / fin de ligne (`Home` et `End` aussi) |
| `Ctrl+Up` / `Ctrl+Down` | début / fin du fichier |

Le panneau de documentation ouvert par `^D` suit le curseur : placé sur `kn`, il
explique `kn`. Quand le panneau est fermé, la barre d'état signale
`^D documents '<mot>'` dès que le mot sous le curseur est documenté.

Il n'y a pas de touche « copier » : `^K` coupe, un premier `^U` remet le texte en
place, un second `^U` ailleurs en fait une copie. `Shift+PgUp` / `Shift+PgDn`
peuvent être interceptés par le terminal (défilement de son historique).

Coloration
----------

* mot-clé connu : bleu ;
* expression `$ ... $` (voir `define`) : magenta, espaces compris, de sorte qu'un
  `/` dans une expression n'est pas pris pour un commentaire ;
* un mot commençant par `#`, `/` ou `!` : vert jusqu'à la fin de la ligne, comme
  le fait `Lhyphen::loadCONF()`.

D'où viennent les mots-clés et la documentation
-----------------------------------------------

Tout ce que l'éditeur sait du langage de l-hyphen est dans `lhyphen.lang` :
mots-clés, documentation et snippets. Documenter un nouveau mot-clé ne demande
donc aucune recompilation. **Quand un mot-clé est ajouté dans `loadCONF`, l'ajouter
aussi ici** (en plus de `DOCUMENTATION.md` et `cheatsheets/`).

Le fichier est cherché, dans l'ordre :

1. `$LHYPHEN_LANG`
2. le répertoire courant
3. à côté de l'exécutable, puis dans ses deux répertoires parents
4. `~/.lhyphen`

Chacun de ces répertoires est essayé directement et via un sous-répertoire
`lhyphen-edit/`, si bien que le `lhedit` copié à la racine de l-hyphen, lancé
depuis la racine (`./lhedit`) ou depuis un répertoire d'exemple (`../../lhedit`),
trouve son fichier. Pour
installer l'éditeur ailleurs, copier `lhyphen.lang` à côté du binaire ou dans
`~/.lhyphen/`.

Son format est décrit dans son en-tête. En bref :

~~~
[keyword] dt
doc: dt [(double) valeur]

  Pas de temps.

[snippet] Temps
dt _valeur_
nstep _valeur_
~~~

Notes
-----

* Le fichier est tenu en mémoire comme un vecteur de lignes et seules les lignes
  visibles sont dessinées : un gros fichier `conf` s'ouvre instantanément.
* Le texte est en UTF-8 (les commentaires accentués sont fréquents) : le curseur,
  l'effacement et l'affichage avancent d'un caractère à la fois. Les caractères
  « larges » (CJK, emoji) occupent en revanche deux colonnes à l'écran et
  décaleraient le curseur.
* `^C` ne quitte pas (les signaux sont désactivés en mode brut) ; utiliser `^X`.

Organisation des sources
------------------------

| Fichier | Rôle |
|---|---|
| `main.cpp` | l'éditeur : mise en page, touches, coloration |
| `text_buffer.hpp` | lignes, curseur, édition, annuler/refaire, recherche |
| `terminal.hpp` | mode brut, décodage des touches, sortie d'une image écran |
| `utf8.hpp` | frontières de caractères et colonnes à l'écran |
| `lhyphen_lang.hpp` | lecture de `lhyphen.lang` |
| `lhyphen.lang` | mots-clés, documentation et snippets |


# photon_FCEL 
(English bellow)
Programme Python/Jetphox pour afficher les section efficaces, les RpA et plus de photons directs et des photons de fragmentation dans les processus $2\rightarrow2$.

Language utilisé: Python, version 3.12

Modules Python utilisés: 
 * numpy
 * matplolib
 * scipy
 * lhapdf (voir https://www.lhapdf.org/ pour l'installation, ou alors utiliser homebrew: https://brew.sh/)
 * uproot

Version de Jetphox utilisé: v1.3.1_4 (https://lapth.cnrs.fr/PHOX_FAMILY/jetphox.html)

Autres logiciels utiles pour la modification:
* ROOT (version 6.36.04 : https://root.cern)

Le fichier *template.pl* permet de calculer automatiquement les diagrammes de fragmentation (symétriques) LO séparément (les 1,3,5,7,8,9,15,16). Il suffit de modifier le bin en $p_\bot$ dans *parameter_template.indat* et la coupure en $p_\bot$ dans *param_histo_template.indat* pour faire votre propre set et de le placer dans un copy_pawres avec la bonne teminaison (e.g. copy_pawres_5 qui représente le dossier des .root pour $p_\bot=5$ GeV).

Avant de lancer le fichier *LO_frag_ratios.C* avec root, il faut calculer les 4 diagrammes non-symétriques 13a/13b/14a/14b. La démarche à suivre est expliquée dans le fichier *Rpa_hadron_frag.py* aux lignes 111-119. Dans ce cas, pour chaque diagramme, il faut rentrer manuellement (par exemple pour 13a, en n'oubliant pas de modifier *parameter.indat* et *param_histo.indat* en concéquence, avec la nomenclature "LO_diag"+n_processus pour le nom du root et de l'éxecutable):

*perl start.pl*

*./runLO_diag13a*

Dans votre dossier avec tous vos .root, vous lancez la commande: 

*root LO_frag_ratios.C*

Vous pouvez maintenant utiliser les fonctions de plots.

# English version

Python/Jetphox program to display cross sections, RpA, and more for direct photons and fragmentation photons in $2\rightarrow2$ processes.

Language used: Python, v3.12

Python modules used: 
 * numpy
 * matplolib
 * scipy
 * lhapdf (see https://www.lhapdf.org/ for installation, or use homebrew: https://brew.sh/ )
 * uproot

Jetphox version used: v1.3.1_4 ( https://lapth.cnrs.fr/PHOX_FAMILY/jetphox.html )

Other useful software for modification:
* ROOT (version 6.36.04 : https://root.cern)

The template.pl file allows you to automatically compute the LO symmetric fragmentation diagrams separately (1,3,5,7,8,9,15,16).
You just need to modify the $p_\bot$ bin in parameter_template.indat and the $p_\bot$ cut in param_histo_template.indat to create your own set, and place it in a copy_pawres folder with the proper suffix (e.g. copy_pawres_5, which represents the folder containing the .root files for $p_\bot = 5$ GeV).

Before running LO_frag_ratios.C with ROOT, you need to compute the 4 non-symmetric diagrams 13a/13b/14a/14b.
The procedure is explained in Rpa_hadron_frag.py between lines 111–119.
In this case, for each diagram, you must run manually (for example for 13a, remembering to modify parameter.indat and param_histo.indat accordingly, using the naming convention “LO_diag” + process_number for the ROOT file and executable):

*perl start.pl*

*./runLO_diag13a*

In the folder containing all your .root files, run the command: 

*root LO_frag_ratios.C* 

You can now use the plotting functions.

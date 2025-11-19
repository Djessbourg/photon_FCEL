
# photon_FCEL
Programme Python/Jetphox pour afficher les section efficaces, les RpA et plus de photons directs et des photons de fragmentation dans les processus 2->2.

Language utilisé: Python, version 3.12

Modules Python utilisés: 
 * numpy
 * matplolib
 * scipy
 * lhapdf (voir https://www.lhapdf.org/ pour l'installation, ou alors utiliser homebrew: https://brew.sh/ )
 * uproot

Version de Jetphox utilisé: v1.3.1_4 ( https://lapth.cnrs.fr/PHOX_FAMILY/jetphox.html )

Autres logiciels utiles pour la modification:
* ROOT ( version 6.36.04 : https://root.cern)

Le fichier *template.pl* permet de calculer automatiquement les diagrammes de fragmentation (symétriques) LO séparément (les 1,3,5,7,8,9,15,16). Il suffit de modifier le bin en $p_\bot$ dans *parameter_template.indat* et la coupure en $p_\bot$ dans *param_histo_template.indat* pour faire votre propre set et de le placer dans un copy_pawres avec la bonne teminaison (e.g. copy_pawres_5 qui représente le dossier des .root pour $p_\bot=5$ GeV).
Avant de lancer le fichier *LO_frag_ratios.C* avec root, il faut calculer les 4 diagrammes non-symétriques 13a/13b/14a/14b. La démarche à suivre est expliquée dans le fichier *Rpa_hadron_frag.py* aux lignes 111-119. Dans ce cas, pour chaque diagramme, il faut rentrer manuellement (par exemple pour 13a, en n'oubliant pas de modifier *parameter.indat* et *param_histo.indat* en concéquence, avec la nomenclature "LO_diag"+n_processus pour le nom du root et de l'éxecutable):

*perl start.pl*

*./runLO_diag13a*

Dans votre dossier avec tous vos .root, vous lancez la commande: *root LO_frag_ratios.C* et vous pouvez maintenant utiliser les fonctions de plots.

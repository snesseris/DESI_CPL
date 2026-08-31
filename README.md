# To CPL, or not to CPL? What we have not learned about the dark energy equation of state
Here we present the MCMC codes and chains used in our analysis for [arXiv: 2503.22529](https://arxiv.org/abs/2503.22529), which has been published in [PRD Letters](https://journals.aps.org/prd/abstract/10.1103/q2gb-n61g).

The repo contains the following:
1) In /data there are the Pantheon+ data needed in the analysis,
2) The bigchain_[xxx].zip files contain the chains used in the paper to make the plots, with the plots_CPL_et_al.nb Mathematica notebook (see also MCMC_CPL_Panthon_plus_v1.nb).
3) The chains were created with the script_[xxx].sh wolfram scripts that were submitted via the script_[xxx].sh scripts to a cluster to run.

To run the codes, you need to download everything in the same folder, make sure you have Mathematica (v14+) installed and make sure you edit any paths in the [xxx].wls scripts if you get errors about not finding the files etc.

If you use any of the above, please cite the [paper](https://arxiv.org/abs/2503.22529), here's also the bib file:
@article{Nesseris:2025lke,
    author = "Nesseris, Savvas and Akrami, Yashar and Starkman, Glenn D.",
    title = "{What we have not learned about the dark energy equation of state}",
    eprint = "2503.22529",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    reportNumber = "IFT-UAM/CSIC-25-29",
    doi = "10.1103/q2gb-n61g",
    journal = "Phys. Rev. D",
    volume = "114",
    number = "4",
    pages = "L041305",
    year = "2026"
}

Created by [Savvas Nesseris](https://members.ift.uam-csic.es/savvas.nesseris/), 2026. For any questions drop me a line at savvas.nesseris@csic.es

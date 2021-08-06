function [Diff01,Diff05,Diff1,Diff15,Diff2,Diff100] = MeanCal(DenseDiff01_run1,DenseDiff01_run2,DenseDiff05_run1,DenseDiff05_run2,DenseDiff1_run1,DenseDiff1_run2,DenseDiff15_run1,DenseDiff15_run2,DenseDiff2_run1,DenseDiff2_run2,DenseDiff100_run1,DenseDiff100_run2)


Diff01.uninfected_live = mean([DenseDiff01_run1.uninfected_live;DenseDiff01_run2.uninfected_live]);
Diff05.uninfected_live = mean([DenseDiff05_run1.uninfected_live;DenseDiff05_run2.uninfected_live]);
Diff1.uninfected_live = mean([DenseDiff1_run1.uninfected_live;DenseDiff1_run2.uninfected_live]);
Diff15.uninfected_live = mean([DenseDiff15_run1.uninfected_live;DenseDiff15_run2.uninfected_live]);
Diff2.uninfected_live = mean([DenseDiff2_run1.uninfected_live;DenseDiff2_run2.uninfected_live]);
Diff100.uninfected_live = mean([DenseDiff100_run1.uninfected_live;DenseDiff100_run2.uninfected_live]);

Diff01.dead = mean([DenseDiff01_run1.dead;DenseDiff01_run2.dead]);
Diff05.dead = mean([DenseDiff05_run1.dead;DenseDiff05_run2.dead]);
Diff1.dead = mean([DenseDiff1_run1.dead;DenseDiff1_run2.dead]);
Diff15.dead = mean([DenseDiff15_run1.dead;DenseDiff15_run2.dead]);
Diff2.dead = mean([DenseDiff2_run1.dead;DenseDiff2_run2.dead]);
Diff100.dead = mean([DenseDiff100_run1.dead;DenseDiff100_run2.dead]);


Diff01.extracellular_virus = mean([DenseDiff01_run1.extracellular_virus;DenseDiff01_run2.extracellular_virus]);
Diff05.extracellular_virus = mean([DenseDiff05_run1.extracellular_virus;DenseDiff05_run2.extracellular_virus]);
Diff1.extracellular_virus = mean([DenseDiff1_run1.extracellular_virus;DenseDiff1_run2.extracellular_virus]);
Diff15.extracellular_virus = mean([DenseDiff15_run1.extracellular_virus;DenseDiff15_run2.extracellular_virus]);
Diff2.extracellular_virus = mean([DenseDiff2_run1.extracellular_virus;DenseDiff2_run2.extracellular_virus]);
Diff100.extracellular_virus = mean([DenseDiff100_run1.extracellular_virus;DenseDiff100_run2.extracellular_virus]);

end
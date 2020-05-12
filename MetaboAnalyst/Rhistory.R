mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<-c("HMDB31923","HMDB00872","HMDB00161","HMDB11532","HMDB00714","HMDB00152","HMDB00671","HMDB02203","HMDB00182","HMDB00222","HMDB02320","HMDB0061116","HMDB00289","HMDB11621","HMDB03073","HMDB13324","HMDB11753","HMDB00094","HMDB0011635","HMDB03374","HMDB13164","HMDB00430","HMDB00623","HMDB02329","HMDB00148","HMDB00782","HMDB0000619","HMDB00562")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-PerformDetailMatch(mSet, "HMDB0061116");
mSet<-GetCandidateList(mSet);
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateOraScore(mSet, "dgr", "fisher")
mSet<-PlotPathSummary(mSet, "path_view_0_", "png", 72, width=NA)
mSet<-PlotKEGGPath(mSet, "D-Arginine and D-ornithine metabolism",528, 480, "png", NULL)
mSet<-RerenderMetPAGraph(mSet, "zoom1589292144657.png",528.0, 480.0, 100.0)
mSet<-PlotKEGGPath(mSet, "D-Arginine and D-ornithine metabolism",528, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "D-Glutamine and D-glutamate metabolism",528, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "Nitrogen metabolism",528, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "D-Arginine and D-ornithine metabolism",528, 480, "png", NULL)

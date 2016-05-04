#require(rPython)
#require(mefa)
#require(ape)

python.load("/afd-python/SourceCode.py")
python.load("/afd-python/suffix_tree.py")
python.exec("import sys")
python.exec("sys.path.append(\"/afd-python/\")")

do_FFP <- function(fasta, k=5, ffppath, dna_names){
  #dist<-python.call("computeffpvectoranddistances",fasta,k,ffppath,temppath)
  sh <- paste0('sudo ', ffppath, 'ffpry -l ', k, ' -d -m ', fasta, ' | ', ffppath, 'ffpcol -d | ', ffppath, 'ffprwn | ', ffppath, 'ffpjsd')
  print(sh)
  distvec <- system(sh, intern = TRUE)
  #dist2<-as.matrix(vec2dist(dist[[2]],length(dist[[1]])))
  #colnames(dist2)<-rownames(dist2)<-dist[[1]]
  #return(dist2)
  distvec2 <- do.call(rbind, lapply(distvec, function(x) as.numeric(strsplit(x, " ")[[1]])))
  rownames(distvec2) <- dna_names
  colnames(distvec2) <- dna_names
  distvec2
}

doCVV<-function(path,temppath,rep,k=5,ffppath,outpath){
  dist<-python.call("computeCVvectorDistances",paste(path,rep,"_All.fa",sep=""),k,ffppath,temppath)
  dist2<-as.matrix(vec2dist(dist[[2]],length(dist[[1]])))
  colnames(dist2)<-rownames(dist2)<-dist[[1]]
  write.csv(dist2,file=paste(outpath,rep,"_K_",k,"_CVV_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_K_",k,"_CVV_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}

doRTD<-function(path, k=5){
  file<-paste(path,rep,"_All.fa",sep="")
  python.exec(paste("(names,vectormatrix)=computeRTDfromfastafile(\"",file,"\",",k,")",sep=""))
  python.exec("import scipy.spatial.distance")
  python.exec("distancematrix=scipy.spatial.distance.pdist(vectormatrix,'euclidean')") #### methoditem include: 'euclidean','braycurtis','canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine','minkowski','seuclidean' or 'sqeuclidean'
  python.exec("distancematrix2=distancematrix.tolist()")
  dist<-python.get("distancematrix2")
  nn<-python.get("names")
  dist2<-vec2dist(dist,length(nn),nn)
  return(dist2)
}

doFCGR<-function(path,rep,k=5,outpath){
  file<-paste(path,rep,"_All.fa",sep="")
  python.exec(paste("(names,vectormatrix)=computeFCGRfromfastafile(\"",file,"\",",k,")",sep=""))
  python.exec("import scipy.spatial.distance")
  python.exec("distancematrix=scipy.spatial.distance.pdist(vectormatrix,'euclidean')") #### methoditem include: 'euclidean','braycurtis','canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine','minkowski','seuclidean' or 'sqeuclidean'
  python.exec("distancematrix2=distancematrix.tolist()")
  dist<-python.get("distancematrix2")
  nn<-python.get("names")
  dist2<-vec2dist(dist,length(nn),nn)
  write.csv(as.matrix(dist2),file=paste(outpath,rep,"_K_",k,"_FCGR_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_K_",k,"_FCGR_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}

doBBC<-function(path,rep,k=5,outpath){
  file<-paste(path,rep,"_All.fa",sep="")
  python.exec(paste("(names,vectormatrix)=computeBBCfromfastafile(\"",file,"\",",k,")",sep=""))
  python.exec("import scipy.spatial.distance")
  python.exec("distancematrix=scipy.spatial.distance.pdist(vectormatrix,'euclidean')") #### methoditem include: 'euclidean','braycurtis','canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine','minkowski','seuclidean' or 'sqeuclidean'
  python.exec("distancematrix2=distancematrix.tolist()")
  dist<-python.get("distancematrix2")
  nn<-python.get("names")
  dist2<-vec2dist(dist,length(nn),nn)
  write.csv(as.matrix(dist2),file=paste(outpath,rep,"_K_",k,"_BBC_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_K_",k,"_BBC_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}

doICPIC<-function(path,rep,k=5,outpath){
  file<-paste(path,rep,"_All.fa",sep="")
  python.exec(paste("(names,vectormatrix)=computeICPICfromfastafile(\"",file,"\",",k,")",sep=""))
  python.exec("import scipy.spatial.distance")
  python.exec("distancematrix=scipy.spatial.distance.pdist(vectormatrix,'euclidean')") #### methoditem include: 'euclidean','braycurtis','canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine','minkowski','seuclidean' or 'sqeuclidean'
  python.exec("distancematrix2=distancematrix.tolist()")
  dist<-python.get("distancematrix2")
  nn<-python.get("names")
  dist2<-vec2dist(dist,length(nn),nn)
  write.csv(as.matrix(dist2),file=paste(outpath,rep,"_K_",k,"_ICPIC_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_K_",k,"_ICPIC_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}

do2DSV<-function(path,rep,outpath){
  file<-paste(path,rep,"_All.fa",sep="")
  python.exec(paste("(names,vectormatrix)=compute2DGraphfromfastafile(\"",file,"\")",sep=""))
  python.exec("import scipy.spatial.distance")
  python.exec("distancematrix=scipy.spatial.distance.pdist(vectormatrix,'euclidean')") #### methoditem include: 'euclidean','braycurtis','canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine','minkowski','seuclidean' or 'sqeuclidean'
  python.exec("distancematrix2=distancematrix.tolist()")
  dist<-python.get("distancematrix2")
  nn<-python.get("names")
  dist2<-vec2dist(dist,length(nn),nn)
  write.csv(as.matrix(dist2),file=paste(outpath,rep,"_2DSV_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_2DSV_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}

do2DMV<-function(path,rep,k=5,outpath){
  file<-paste(path,rep,"_All.fa",sep="")
  python.exec(paste("(names,vectormatrix)=compute2DGraphMfromfastafile(\"",file,"\",",k,")",sep=""))
  python.exec("import scipy.spatial.distance")
  python.exec("distancematrix=scipy.spatial.distance.pdist(vectormatrix,'euclidean')") #### methoditem include: 'euclidean','braycurtis','canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine','minkowski','seuclidean' or 'sqeuclidean'
  python.exec("distancematrix2=distancematrix.tolist()")
  dist<-python.get("distancematrix2")
  nn<-python.get("names")
  dist2<-vec2dist(dist,length(nn),nn)
  write.csv(as.matrix(dist2),file=paste(outpath,rep,"_K_",k,"_2DMV_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_K_",k,"_2DMV_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}

do2DNV<-function(path,rep,outpath){
  file<-paste(path,rep,"_All.fa",sep="")
  python.exec(paste("(names,vectormatrix)=computeNatureVecfromfastafile(\"",file,"\")",sep=""))
  python.exec("import scipy.spatial.distance")
  python.exec("distancematrix=scipy.spatial.distance.pdist(vectormatrix,'euclidean')") #### methoditem include: 'euclidean','braycurtis','canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine','minkowski','seuclidean' or 'sqeuclidean'
  python.exec("distancematrix2=distancematrix.tolist()")
  dist<-python.get("distancematrix2")
  nn<-python.get("names")
  dist2<-vec2dist(dist,length(nn),nn)
  write.csv(as.matrix(dist2),file=paste(outpath,rep,"_2DNV_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_2DNV_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}

doACS<-function(path,rep,outpath){
  dist<-python.call("computeACSfromfastafile",paste(path,rep,"_All.fa",sep=""))
  dist2<-as.matrix(vec2dist(dist[[2]],length(dist[[1]])))
  colnames(dist2)<-rownames(dist2)<-dist[[1]]
  write.csv(dist2,file=paste(outpath,rep,"_ACS_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_ACS_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}

doKR<-function(path,rep,outpath,temppath,krpath){
  dist<-python.call("computeKrdistances",paste(path,rep,"_All.fa",sep=""),krpath,temppath)
  dist2<-as.matrix(vec2dist(dist[[2]],length(dist[[1]])))
  colnames(dist2)<-rownames(dist2)<-dist[[1]]
  write.csv(dist2,file=paste(outpath,rep,"_KR_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_KR_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}

doLZ<-function(path,rep,outpath,temppath,lzpath){
  dist<-python.call("computeLZdistances",paste(path,rep,"_All.fa",sep=""),lzpath,temppath)
  dist2<-as.matrix(vec2dist(dist[[2]],length(dist[[1]])))
  colnames(dist2)<-rownames(dist2)<-dist[[1]]
  write.csv(dist2,file=paste(outpath,rep,"_LZ_Dist.csv",sep=""))
  tree<-nj(dist2)
  write.tree(tree,file=paste(outpath,rep,"_LZ_NJ_Tree.tre",sep=""))
  out<-list(Tree=tree,Distance=as.dist(dist2))
  return(out)
}
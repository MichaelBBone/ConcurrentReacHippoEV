#### Item Specific Reactivation ####
# assumes all required files are in the working directory
#
# Neural activity predictions are generated/saved in section "Average Encoding Activity Per Image"
# Predictions are used to generate reactivation rank values in section "Reactivation Rank Measure"
#
# ROIs divided into 'chunks' (see "Get Reactivation Results by ROI 'chunk'" section). Only data for 'LowOcc' and 'Hippo' are included.
#
# To generate/save the reactivation values simply run the whole script. If you just want to get data for one ROI chunk, modify 'useChunks'.
#
# Activity predictions (avActSub) and are saved in:
# "avActRoi<ROI chunk>.RData". The data is a list (names = subject) of 2-D arrays (row = image (alpha order), column = vertex)
#
# Reactivation ranks for recall (rankTriSubPrbROI) and recognition (rankTriSubPrbROI) are saved in 
# "reacRankRoi<ROI chunk>_ItemSpec.RData". The data is a list (names = ROI) of 2-D arrays (row = subject, column = trial/image (alpha order))
#
# hippocampal data split into 5 sections
# ROI names: "17_1" "17_2" "17_3" "17_4" "17_5" "53_1" "53_2" "53_3" "53_4" "53_5"
# 17 is left, 53 is right hippocampus.
# 1 is most posterior, 5 is most anterior


#### Global Variables ####

smooth = T # use spatially smoothed brain data
subjectsEx = c(1002,1006,1007,1008,1009,1010,1011,1012,1014,1016,1019,1020,1022,1025,1026,1027,
               1028,1030,1032,1033,1034,1035,1036,1037,1038) # subject numbers
FeatLvlNames = c('L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12','L13','FC1','FC2','FC3') # feature level names

# get groups of Freesurfer ROIs
# because it takes so long to get the encoding model predictions I did it in 'chunks' of ROIs
load(file='ROINamesDiv.RData')
ROILvlNames = c('LowOc','HighOc','Tmprl','Prtl','Front','Other','Hippo') # names of 'chunks'

# load image names
# nBackImagesAll = all image names
# list (names = subjects) of vectors containing image names:
# nBackImages = encoding image for each trial
# recogImages = cued recall images for each image and time point (sorted alphabetically)
# recogProbes = recognition probe images for each image and time point (sorted alphabetically)
load(file=paste0("imageNames.RData"))



#### Get Reactivation Results by ROI 'chunk' #####
# there are 7 chunks:
# 1 = (LowOc) occipital early visual
# 2 = (HighOc) occipital late visual
# 3 = (Tmprl) temporal
# 4 = (Prtl) parietal
# 5 = (Front) frontal
# 6 = (Other) every other neocortical ROI
# 7 = (Hippo) hippocampus
useChunks = c(1,7)

ROILvl = 7
for (ROILvl in useChunks) {
  print(ROILvlNames[ROILvl])
  tROINames = ROINamesDiv[[ROILvl]]
  
  # load ROI name of each vertex (ROINodeNames)
  # list (names = subjects) of vectors containing roi names
  load(file=paste0("roiNodeNames",ROILvlNames[ROILvl],".RData"))
  
  # load movie neural data (movie2VntrDatSep, movie3VntrDatSep)
  # list (names = subjects) of matrices (row = time point, column = brain vertex)
  # data z-scored for each vertex
  load(file=paste0("movie",ROILvlNames[ROILvl],".RData"))
  
  # load encoding neural data (nBackVntrDat)
  # list (names = subjects) of matrices (row = trial, column = brain vertex)
  # data z-scored for each vertex
  load(file=paste0("nback",ROILvlNames[ROILvl],".RData"))
  
  # load recall (recogVis) and recognition (recogPrb) neural data
  # list (names = subjects) of matrices (row = image (alphabetically ordered), column = brain vertex)
  load(paste0("recogVis",ROILvlNames[ROILvl],".RData"))
  load(paste0("recogPrb",ROILvlNames[ROILvl],".RData"))
  
  
  #### Average Encoding Activity Per Image ####
  sub = 1002
  avActSub = list()
  rankTriSubPrb = array(NA,c(length(subjectsEx),nrow(recogPrb[[1]])))
  rownames(rankTriSubPrb) = subjectsEx
  rankTriSubVis = rankTriSubPrb
  for (sub in subjectsEx) {
    print(sub)
    # get brain data
    dat = nBackVntrDat[[toString(sub)]]
    datPrb = scale(recogPrb[[toString(sub)]])
    datVis = scale(recogVis[[toString(sub)]])
    
    # get encoding image data
    datNames = nBackImages[[toString(sub)]]
    tempIX = sapply(datNames,function (x) {which(nBackImagesAll==x)})
    imOrder = ceiling(tempIX/2)
    
    # get average activity per image
    avAct = array(NA,c(90,ncol(dat)))
    im = 1
    for (im in 1:90) {
      avAct[im,] = colMeans(dat[which(imOrder==im),])
    }
    avActSub[[toString(sub)]] = avAct
    
    # print reactivation for region (during probe and visualization--mostly interested in visualization)
    voxIX = 1:ncol(dat)
    actCor = cor(t(avActSub[[toString(sub)]][,voxIX]),t(datPrb[,voxIX]))
    for (trial in 1:nrow(datPrb)) {
      rankTriSubPrb[toString(sub),trial] = which(sort.int(actCor[trial,],decreasing=T,index.return=T)$ix==trial)
    }
    print(45.5-mean(rankTriSubPrb[toString(sub),]))
    
    actCor = cor(t(avActSub[[toString(sub)]][,voxIX]),t(datVis[,voxIX]))
    for (trial in 1:nrow(datVis)) {
      rankTriSubVis[toString(sub),trial] = which(sort.int(actCor[trial,],decreasing=T,index.return=T)$ix==trial)
    }
    print(45.5-mean(rankTriSubVis[toString(sub),]))
  }
  save(avActSub,file=paste0("avActRoi",ROILvlNames[ROILvl],".RData"))
  
  
  #### Reactivation Rank Measure ####
  # generates reactivation rank as described in the paper
  
  # get FreeSurfer ROI names in the current larger ROI
  uniqueROIs = sort(unique(ROINodeNames[[1]]))
  
  # chance rank
  rankChance = mean(1:90)
  
  # load activity prediction
  datVars = load(file=paste0("avActRoi",ROILvlNames[ROILvl],".RData"))
  
  # initialize data structures
  rankTriSubPrbROI = list()
  rankTriSubVisROI = list()
  corDifPrbVsOthrSubVisROI = list()
  TValROI = rep(NA,length(uniqueROIs))
  names(TValROI) = uniqueROIs
  
  ## get reactivation ranks for each FreeSurfer ROI
  roi = uniqueROIs[4]
  sub = 1002
  countROI = 0
  for (roi in uniqueROIs) {
    countROI = countROI + 1
    print(roi)
    
    # initialize data structures
    rankTriSubPrb = array(NA,c(length(subjectsEx),nrow(recogPrb[[1]])))
    rownames(rankTriSubPrb) = subjectsEx
    rankTriSubVis = rankTriSubPrb
    
    ## get reactivation ranks for each subject
    sub = subjectsEx[1]
    count = 0
    for (sub in subjectsEx) {
      count = count + 1
      print(sub)
      
      # get column vector index for current ROI
      ROIIx = which(ROINodeNames[[toString(sub)]]==roi)
      
      # get recall (vis) and recognition (prb) brain data for subject
      datVis = scale(recogVis[[toString(sub)]])
      datPrb = scale(recogPrb[[toString(sub)]])
      
      # get encoding model predictions
      predDatVis = eval(parse(text=paste0(datVars,"[[toString(sub)]]")))
      
      # get number of vertices and images
      voxN = ncol(datPrb)
      imN = nrow(datPrb)
      
      # limit vertices to those in the target ROI
      voxIX = ROIIx
      
      # generate/print trial-average reactivation (rank) results during recall for the current ROI:feature-level:subject
      predCor = cor(t(datVis[,voxIX]),t(predDatVis[,voxIX]))
      for (trial in 1:nrow(datPrb)) {
        rankTriSubVis[toString(sub),trial] = which(sort.int(predCor[trial,],decreasing=T,index.return=T)$ix==trial)
      }
      
      # generate/print trial-average reactivation (rank) results during recognition for the current ROI:feature-level:subject
      predCor = cor(t(datPrb[,voxIX]),t(predDatVis[,voxIX]))
      for (trial in 1:nrow(datPrb)) {
        rankTriSubPrb[toString(sub),trial] = which(sort.int(predCor[trial,],decreasing=T,index.return=T)$ix==trial)
      }
    }
    rankTriSubVisROI[[toString(roi)]] = rankChance - rankTriSubVis
    rankTriSubPrbROI[[toString(roi)]] = rankChance - rankTriSubPrb
  }
  
  # save rank data
  save(rankTriSubPrbROI,rankTriSubVisROI,
       file=paste0("reacRankRoi",ROILvlNames[ROILvl],"_ItemSpec.RData"))
}

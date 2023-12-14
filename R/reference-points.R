#### Reference values for various things

##### COORDONATE SYSTEMS #####

cn.gc <- function(){
  #' Granite coordinates
  #' Granite coords are defined by (Si,Al,Ca,NK = Na+K,FM = Fe+Mg)

  return(c("Si","Al","Ca","NaK","FM"))
}

##### Some useful objects #####
# A = CNK plane
acnk <- function(coordSys="gc"){
  #' Returns the values for the A/CNK = 1 plane in the appropriate coord sys
  #' @param coordsys coordinate system to use.
  #'    gc: granite coordinates (Si,Al,Ca,NK = Na+K,FM = Fe+Mg)


  if(coordSys=="gc"){
    ee <- c(3,1,0,1,0,	#fsp
            2,2,1,0,0,	#an
            0,0,0,0,1,  # top
            3,1,0,1,0) # back to fsp

    nn <-c("fsp","an","FM","fsp")

    prepareData(ee,cn.gc(),nn)

  }
}


# Ideal minerals
idmins <- function(coordSys="gc"){
  #' Returns the values for some ideal mineral compositions in the appropriate coord sys
  #' @param coordsys coordinate system to use.
  #'    gc: granite coordinates (Si,Al,Ca,NK = Na+K,FM = Fe+Mg)


  if(coordSys=="gc"){
    im<-c(1,0,0,0,0,	#q
          3,1,0,1,0,	#fsp
          2,2,1,0,0,	#an
          2.5,1.5,0.5,0.5,0,	#an50
          3,3,2,0,0,	#cz
          3,2,2,0,1,	#Ep
          1,1,0,0,0,	#sill
          2,0,0,0,2,	#opx
          2,0,1,0,1,	#cpx
          1,0,0,0,2,	#olv
          3,2,3,0,0,	#grs
          3,2,0,0,3,	#Gt
          5,4,0,1,2,	#NaCrd
          5,4,0,0,2,	#Crd
          6,2,0,2,6,	#bio
          6,6,0,2,0,	#ms
          7,2,2,0,4,	#MgHbl
          7,1,2,1,5,	#Edn
          6,3,2,1,4)	#Pgs

    nn <-c("q","fsp","an","an50","cz","Ep","sill","opx","cpx","olv","grs-Gt","Gt","NaCrd","Crd","bio","ms","MgHbl","Edn","Pgs")

    prepareData(im,cn.gc(),nn,Mineral)

  }
}

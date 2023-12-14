library(dplyr)
library(tidyr)
library(readxl)

# library(tetraCoords)

####Read some data

read_excel("U:\\Recherche\\@Databases\\all_experiments.xlsx") %>%
  replace_na(list(FeO=0,Fe2O3=0)) %>%
  mutate(FeOt = FeO + Fe2O3/1.111) %>%
  # We first calculate millications
  mutate(Si.m = SiO2 / 60.078 * 1000,
         Ti.m = TiO2 / 79.865 * 1000,
         Al.m = Al2O3 / 101.961 * 2 * 1000,
         Fe.m = FeOt / 71.839 * 1000,
         Mn.m = MnO / 70.937 * 1000,
         Mg.m = MgO / 40.299 * 1000,
         Ca.m = CaO / 56.077 * 1000,
         Na.m = Na2O / 61.979 * 2 * 1000,
         K.m = K2O / 94.195 * 2 * 1000,
         Al = Al.m,
         Ca = Ca.m,
         NaK = Na.m + K.m,
         FM = Fe.m + Mg.m) %>%
  {.} -> dd

# Prepare it for plotting in tetrahedron

dd %>%
  mutate(tetraCoords( Al,Ca,NaK,FM ) ) %>%
  {.} -> dd.gcn

# A/CNK plane
ACNKplane <- acnk() %>% mutate(tetraCoords( Al,Ca,NaK,FM ) )

# Ideal minerals
idMins <- idmins() %>% mutate(tetraCoords( Al,Ca,NaK,FM ) )

# Some plotting
p2 <- drawggTetraPlot(apicesNames=c("Al","Ca","NaK","FM"),show=F) +
  geom_sphere_3d(data=dd.gcn,aes(x=x,y=y,z=z,colour=Colour))+
  scale_colour_identity() +
  geom_path_3d(data=ACNKplane, aes(x=x,y=y,z=z))+
  geom_sphere_3d(data=idMins,aes(x=x,y=y,z=z),color="red")+
  geom_text_z(data=idMins,aes(x=x,y=y,z=z,label=Mineral))


# 3D interactive plot
devoutrgl::rgldev(zscale=3.5)
p2

# For export to webGL, we must capture the graph device to a regular RGL device
# Careful, webGL becomes ~unusable with > 500 points or so
open3d()
p2
aspect3d(1,1,sqrt(6)/3)

s<-scene3d()
rglwidget(s,elementId="tetraWidget")


### Projection from biotite

pb <- c(3,0,2,0, #ms1
        1,0,1,0, #fsp
        1,1,0,0, #CaAl
        1,0,1,3) #bio

cn <-c("Al","Ca","NaK","FM")
nn <- c("xms1","xfsp","xCaAl","xbio")

projPoles <- prepareData(pb,cn,nn,"endMembers")  %>% mutate(tetraCoords( Al,Ca,NaK,FM ) )
tm <-        dataTib2Mat(projPoles, cn,"endMembers")
mydataProjected<-coordMap(dd,tm) %>% mutate(tetraCoords( xms1,xfsp,xCaAl,xbio ) )


## 3D (regular)
p3 <- drawggTetraPlot(apicesNames=c("Al","Ca","NaK","FM"),show=F) +
  geom_sphere_3d(data=dd.gcn,aes(x=x,y=y,z=z,colour=Colour))+
  scale_colour_identity() +
  geom_path_3d(data=ACNKplane, aes(x=x,y=y,z=z))+
  geom_sphere_3d(data=idMins,aes(x=x,y=y,z=z),color="black",size=3)+
  geom_text_z(data=idMins,aes(x=x,y=y,z=z,label=Mineral))+
  geom_sphere_3d(data=projPoles,aes(x=x,y=y,z=z),color="blue",size=10)


# 3D interactive plot
devoutrgl::rgldev(zscale=3.5)
p3

#3D (expanded as in projBiot)

# ACNK plane, clipped
ee <- c(3,1,0,1,0,	#fsp
        2.5,1.5,.5,.5,0,	#an50
        6,2,0,2,6,  # bio
        3,1,0,1,0) # back to fsp

nnm <-c("fsp","an","FM","fsp")
aa <- prepareData(ee,cn.gc(),nnm) %>% coordMap(tm) %>% mutate(tetraCoords( xms1,xfsp,xCaAl,xbio ) )

# project minerals
idMins2 <- idmins() %>%
  coordMap(tm) %>%
  mutate(tetraCoords( xms1,xfsp,xCaAl,xbio ) ) %>%
  filter( abs(x)<1&abs(y)<1&abs(z)<1&!is.na(x))

p4 <- drawggTetraPlot(apicesNames=nn,show=F) +
  geom_sphere_3d(data= (mydataProjected %>% filter(abs(x)<1&abs(y)<1&abs(z)<1)) ,aes(x=x,y=y,z=z,colour=Colour))+
  scale_colour_identity() +
  geom_path_3d(data=aa, aes(x=x,y=y,z=z))+
  geom_sphere_3d(data=idMins2,aes(x=x,y=y,z=z),color="black",size=3)+
  geom_text_z(data=idMins2,aes(x=x,y=y,z=z,label=Mineral))

# 3D interactive plot
devoutrgl::rgldev(zscale=1.5)
p4



# flat Biotite proj diagram

# Subset...
# foo<- dd %>% head(5) %>% select(1,"SiO2",cn,"ms1","fsp","CaAl","bio")



# Ploting coordinates
mydataProjected %<>%
mutate(ss = ms1 + fsp + CaAl,
       aa = CaAl / ss,
       bb = ms1 / ss,
       cc = fsp / ss,
       xcoord = sqrt(3)/2*cc,
       ycoord = (bb-aa)/2 )


mydataProjected %>%
  # Basic plot definition
  ggplot(aes(x=xcoord,y=ycoord))+
  # We plot points
  # Color and Symbol are taken from the file (i.e. the columns must exist, or else !)
  geom_point(aes(color=`Src_type(Jensen)`))+
  # We use the col and pch values found in the file, no mapping is done
  #scale_color_identity()+
  # Draw the decorations
  annotate("path",x=c(-2,1),y=c(0,0),color="black")+ # Horiz line
  annotate("path",x=c(0,0,sqrt(3)/2,0),y=c(-1/2,1/2,0,-1/2),color="black")+ # Triangle
  annotate("path",x=c(-2,sqrt(3)/2),y=c(1/6*(1+2/(sqrt(3)/2)),0),col="grey",lty="dashed" )+ # Feldspar line
  # Text
  annotate("text",x=0,y=-0.55,label="Ca+Al",adj=0.5)+
  annotate("text",x=0,y=0.55,label="3 Al+2(Na+K)",adj=0.5)+
  annotate("text",x=sqrt(3)/2+.03,y=0+0.03,label="Al+(Na+K)",adj=0)+
  # We fix the plot boundaries, and critically set aspect ratio to 1
  coord_fixed(ratio=1,xlim=c(-2,1),ylim=c(-1/2-0.1,1/2+0.1))



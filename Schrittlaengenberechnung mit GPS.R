work_dir <- "C:/Users/joma9/OneDrive/Desktop/RStudio_Workspace"
setwd(work_dir)

library(readr)
library(pracma)
library(rootSolve)
library(signal)
library(geosphere)

#F¸r find.malibrary()
library(acepack)
library(Hmisc)

#CSV datei einlesen und in dat speichern
dat <- read.csv(file="CSV/Messung 26_08/linear_nichlinear_30gross_50normal_42klein/Accelerometer.csv", sep=",", header=TRUE)[,c(1,2,3,4)]
colnames(dat) <- c("Time", "Acc.X", "Acc.Y", "Acc.Z")
View(dat)

#GPS Werte einlesen
datgps <- read_csv("CSV/Messung 26_08/linear_nichlinear_30gross_50normal_42klein/Location.csv")[,c(1,2,3,7)]
colnames(datgps) <- c("Time", "Lat","Long", "HorAccuracy")


#Wieviele Sekunden dauerte die Messung von dat
t <- (tail(dat$Time, n=1))
#Eintraege insgesamt
anzEintraege <- length(dat$Time)

#dt zwischen Messungen. Entweder 0.01(Iphone, Samsung) oder 0.0025(Xiaomi)
messabstandT <- t/anzEintraege
#Hertz. Entweder 394Hz oder 100Hz
frequenz <- round(1/messabstandT)

#alle Indizes an denen datgps$HorAccuracy kleiner als 13 ist
indexInacc <- which(datgps$HorAccuracy<13)
#datgps wird mit den kleiner als 13 Werten gef¸llt 
datgps <- cbind(datgps[indexInacc,])
rownames(datgps) <- NULL
View(datgps)

#Index finden AB wann dat = datgps ist (datgps wurde gek¸rzt --> HorAccuracy<10)
ergebnisFindMatchesDatAnf <- find.matches(dat$Time, datgps$Time[1], tol=messabstandT)

indexdatInaccAnf <- which(ergebnisFindMatchesDatAnf$matches == 1)
#Index wo dat = (gek¸rztes) datGPS ist
indexdatInaccAnf <- indexdatInaccAnf[1]

#wieviele Sekunden wir am anfang abschneiden
tcutanf <- dat$Time[indexdatInaccAnf]

#Index finden BIS wann dat = datgps ist (datgps wurde gek¸rzt --> HorAccuracy<13)
lastIndexDatgps <- tail(datgps$Time, n=1)
lastIndexDat <- tail(dat$Time, n=1)
if (lastIndexDatgps<lastIndexDat) {
  lastIndex <- lastIndexDatgps
} else {
  #F¸r den Fall, dass datgps l‰nger aufgenommen hat als dat. datgps wird dann auf die L‰nge von dat gek¸rzt
  lastIndex <- lastIndexDat
  tg <- find.matches(datgps$Time, lastIndex, tol=0.9)
  sd <- which(tg$matches == 1)
  indexDatGpsEnd <- sd[1]
  indexbisDatgps <- which(datgps$Time<=datgps$Time[indexDatGpsEnd])
  datgps$Time[indexDatGpsEnd]
  datgps <- cbind(datgps[indexbisDatgps,])
}
ergebnisFindMatchesDatEnd <- find.matches(dat$Time, lastIndex, tol=messabstandT)
indexdatInaccEnd <- which(ergebnisFindMatchesDatEnd$matches == 1)
indexdatInaccEnd <- indexdatInaccEnd[1]


#alles was grˆﬂer als der erste genaue sekundenwert von datgps ist und kleiner als der letzte genaue Wert
#wird jetzt in ein neues dataframe kopiert
indexabInacc <- which(dat$Time>dat$Time[indexdatInaccAnf])
indexbisInacc <- which(dat$Time<dat$Time[indexdatInaccEnd])

#Herausfinden ob am Ende Werte wegen GPS-Ungenauigkeit abgeschnitten wurden.
tcutend <- (lastIndexDat-dat$Time[(tail(indexbisInacc, n=1))])-messabstandT
if (tcutend > 1) {
  tcutend <- tcutend
} else{
  tcutend <- 0
}

dat1 <- cbind(dat[indexbisInacc,])
dat2 <- cbind(dat[indexabInacc,])
#beide dataframes zusammenf¸hren und in dat abspeichern
dat <- merge(dat1, dat2,  by=c("Time"),all=FALSE)
#die doppelten Spalten lˆschen
dat <- subset(dat, select = -c(5, 6,7))
colnames(dat) <- c("Time", "Acc.X", "Acc.Y", "Acc.Z")


#Funktion butterworthfilter aus dem Seminarbegleitendem Buch
lowPassFilter = function(df, col, samplingFrequency, cutoff, nOrder) {
  nyquistFrequency = 0.5 * samplingFrequency 
  W = cutoff/nyquistFrequency 
  butterworthFilter = butter(nOrder,W)
  df[,paste(col,"_lowpass",sep="")] = filtfilt(butterworthFilter,df[,col])
  return(df)
}
#dat-Werte durch den LowPassFilter gl‰tten
dat_lowPassX <- lowPassFilter(dat, "Acc.X", frequenz, 3, 2)
dat_lowPassY <- lowPassFilter(dat, "Acc.Y", frequenz, 3, 2)
dat_lowPassZ <- lowPassFilter(dat, "Acc.Z", frequenz, 3, 2)

#magnitude der Accelerometerdaten ausrechnen
mag <- sqrt((dat_lowPassX$Acc.X_lowpass)^2 + (dat_lowPassY$Acc.Y_lowpass)^2 + (dat_lowPassZ$Acc.Z_lowpass)^2)

#magnitude minus mittelwert von Magnitude
magNoG = mag - mean(mag)
plot(dat$Time, magNoG,  type = "l", col="red", xlab = "Zeit", ylab = "Beschleunigung (Magnitude)")

#Mindesthoehe, die Acc haben muss um als Hochpunkt erkannt zu werden
minPeakHeight = 0.65

#Findet Peaks und Slopes (tiefPunkte)
Peaks <- findpeaks(magNoG, minpeakheight = minPeakHeight)

#Fuer die uniroot Berechnung. Findet den ersten Tiefpunkt  (Spalte 3 = Slopes)
firstSlope <- Peaks[,3][1]

#Fuer Uniroot braucht den Index. id entspricht dem Index des Dataframes
dat$id <- seq.int(nrow(dat))
#wieviele Daten sind im Dataframe dat (f?r uniroot berechnung)
lengthmagNoG <- length(dat$id)
limit <- lengthmagNoG
l <- firstSlope
k <- l+50
x.points <- c()

while(k<=limit){
  #Hiermit werden die Nulldurchg‰nge berechnet
  zerocrossings <- uniroot.all(approxfun(dat$id, magNoG), interval = range(dat$id[l:k]))
  x.points <- c(x.points, zerocrossings)
  l =  k
  k =  k+50
}


limit <- length(x.points)
x <- 1
y <- 3
i <- 1

#Dataframe fuer Schrittlaengen aus Integration
dfPt <- data.frame("Pt")

while(y<=limit) {
  anf <- round(x.points[x])
  end <- round(x.points[y])
  #1. Integration f¸r die Geschwindigkeit
  Vt <- cumtrapz(dat$Time[anf:end], magNoG[anf:end])
  #2. Integration f¸r die Position/Strecke
  Pt <- cumtrapz(dat$Time[anf:end], Vt)
  #letzter Pt Eintrag wird in dfPt eingef¸gt
  lastPt <- (tail(Pt, n=1))
  dfPt[i, ] = c(lastPt)
  
  x <- x+2
  y <- y+2
  i <- i+1
}
#dfPt numerisch, damit wir mean und sd nutzen kˆnnen
dfPtnumeric <- as.numeric(dfPt$X.Pt.)

numSteps <- length(x.points)/2
cat(numSteps, sep = '', " Schritte insgesamt")

StreckeIntegration <- sum(dfPtnumeric)
cat(round(StreckeIntegration, digits = 3),"m", sep = '', " Strecke errechnet mit Integration")

#dflengthswitch erstellen. Dataframe in dem der Zeitpunkt einer Schrittl‰ngenver‰nderung festgehalten wird 
#+ wieviele Schritte jedes Schrittl‰ngenintervall hatte
x <- 1
y <- 2
j <- 3
xplus2 <- 1
yplus2 <- 3
lastIndexdfLengthSwitch <- 0
dflengthswitch <- data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('Index','Time','Steps'))), stringsAsFactors=F)

laststep <- 0
while(y<=numSteps) {
  anf <- round(x.points[xplus2])
  end <- round(x.points[yplus2])
  if((abs(dfPtnumeric[x]-dfPtnumeric[y])>0.07) & (abs(dfPtnumeric[x]-dfPtnumeric[j])>0.07) | (y == numSteps)) {
    Index <- dat$id[end]
    difference <- Index - lastIndexdfLengthSwitch 
    
    if(difference > 500) {
      
      #Findet Peaks und Slopes (tiefPunkte) von anfang bis magNoG[Index]
      PeaksLengthswitch <- findpeaks(magNoG[1:Index], minpeakheight = 0.5)
      #ACHTUNG PeaksLengthswitch[,3] beinhaltet die Slopes. Wir rechnen mit den Werten nicht. Wir z‰hlen nur, wieviele Slopes es gibt
      #also wieviele Schritte es in dem Intervall magNoG[1:Index] gibt
      PeaksLengthswitch[,3]
      stepsat <- (length(PeaksLengthswitch[,3]))
      stepsat
      
      realsteps <- stepsat - laststep
      realsteps
      
      laststep <- stepsat
      laststep
      
      Timeoflengthswitch <- dat$Time[Index]
      print(Timeoflengthswitch)
      #stepsatlengthswitch <- (length(PeaksLengthswitch[,3]))
      newrow = c(Index, Timeoflengthswitch, realsteps)
      names(dflengthswitch) <- c('Index','Time','Steps')
      dflengthswitch <- rbind(dflengthswitch, newrow)
      lastIndexdfLengthSwitch <- Index   
    }
    
  }else{
    #Optional 
    print("keine Schrittl‰ngenver‰nderung bemerkt")
  }
  xplus2 <- xplus2+2
  yplus2 <- yplus2+2
  x <- x+1
  y <- y+1
  j <- j+1
}
#-------------------------------------------------------------------------------------

#AB HIER IST LOCATIONBERECHNUNG
limiter <- length(dflengthswitch$Time)
p <- 1
result <- 0
indexMatchanf <- 1
while(p<=limiter) {
  #findet die Stellen von dflengthswitch$Time in datgps$Time
  ergebnisFindMatches <- find.matches(datgps$Time, dflengthswitch$Time, tol=0.5)
  #index an der Stelle wo datgps und dflengthswitch matchen
  indexMatchend <- which(ergebnisFindMatches$matches == p)
  #berechnet die distanz aufgrund der GPS Daten
  distanzgps <-(distm(c(datgps$Long[indexMatchanf] , datgps$Lat[indexMatchanf]),
                      c(datgps$Long[indexMatchend] , datgps$Lat[indexMatchend]),
                      fun = distHaversine))
  cat(round(distanzgps, digits = 3),sep = ' ',"m (Strecke)\n")
  #write.table(distanzgps,row.names=F, col.names=F)
  cat(round(distanzgps/dflengthswitch$Steps[p], digits = 3), sep = ' ', "m (Schrittl‰nge)\n")
  indexMatchanf <- which(ergebnisFindMatches$matches == p)
  
  #F‹R JONAS&NILS: >0.95 = groﬂ,   0.65-0.95 = normal,     <0.65 = klein 
  if(distanzgps/dflengthswitch$Steps[p]>=0.95) {
    cat(dflengthswitch$Steps[p],sep = ' ',"groﬂe Schritte\n")
  }
  else if(distanzgps/dflengthswitch$Steps[p]<=0.65) {
    cat(dflengthswitch$Steps[p],"m",sep = ' ',"kleine Schritte\n")
  }
  else {
    cat(dflengthswitch$Steps[p],sep = ' ',"normalgroﬂe Schritte\n")
  }
  cat("\n")
  cat("\n")
  p <- p+1
  result <- result + distanzgps
}
cat(tcutanf,sep = '  ',"Sekunden wurden vorne abgeschnitten")
cat(tcutend,sep = '  ',"Sekunden wurden hinten abgeschnitten")
#View(dflengthswitch)

#hier werden die Schrittl‰ngen mittels GPS ermittelt, wenn alle Schritte gleich sind. Die Werte in dem dataframe tempSchrittlaengen zwischengespeichert und am Ende f¸r die Gesamtstreckenl‰nge zusammenaddiert 
tempSchrittlaengen <- NULL;
lim <- nrow(datgps)
x <- 1
y <- 2
result <- 0
while (y<=lim)
{
  z <- (distm(c(datgps$Long[x] , datgps$Lat[x]), c(datgps$Long[y] , datgps$Lat[y]), fun = distHaversine))
  tempSchrittlaengen <- rbind(tempSchrittlaengen, z)
  x=x+1
  y=y+1
  result=result+z
}


#f¸r die Openstreet Map
library(leaflet)
m <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(lng= datgps$Long, lat= datgps$Lat, popup="The birthplace of R")
m  # Print the map


streckeGPS <- result
#Strecke mit GPS errechnet:
cat(streckeGPS, "m", sep = '', " Strecke insgesamt errechnet mit GPS")
#Schrittl‰nge mit GPS:
schrittlaengeGPS <- streckeGPS/numSteps
cat(round(schrittlaengeGPS, digits = 3),"m", sep = '', " durchschnittliche Schrittl‰nge f¸r die ganze Strecke, errechnet mit GPS")
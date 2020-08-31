work_dir <- "C:/Users/joma9/OneDrive/Desktop/RStudio_Workspace"
setwd(work_dir)

library(readr)
library(pracma)
library(rootSolve)
library(signal)
library(geosphere)

#CSV datei einlesen und in dat speichern
dat <- read.csv(file="CSV/Messung_10_08/Brust_100m_50normale_30grosse_40kleine/Accelerometer.csv", sep=",", header=TRUE)[,c(1,2,3,4)]
colnames(dat) <- c("Time", "Acc.X", "Acc.Y", "Acc.Z")
View(dat)

#Wieviele Sekunden dauerte die Messung von dat
t <- (tail(dat$Time, n=1))
#Eintraege insgesamt
anzEintraege <- length(dat$Time)

#dt zwischen Messungen. Entweder 0.01(Iphone, Samsung) oder 0.0025(Xiaomi)
messabstandT <- t/anzEintraege
#Hertz. Entweder 394Hz oder 100Hz
frequenz <- round(1/messabstandT)

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
k <- l+200
x.points <- c()
while(k<=limit){
  #Hiermit werden die Nulldurchg‰nge berechnet
  zerocrossings <- uniroot.all(approxfun(dat$id, magNoG), interval = range(dat$id[l:k]))
  x.points <- c(x.points, zerocrossings)
  l =  k
  k =  k+200
}


#checken ob von x.points[1] bis x.points[2] negative Zahlen sind. Falls ja, dann m¸ssen wir x und y anpassen
rslt <- 0
for (i in 10){
  m <- magNoG[x.points[1]+i]
  rslt <- rslt+m
}
if (rslt<0) {
  x <- 2
  y <- 4
}else{
  x <- 1
  y <- 3
}
i <- 1
limit <- length(x.points)
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
dfPtnumeric <- dfPtnumeric[dfPtnumeric > 0]
dfPtmean <- mean(dfPtnumeric)
dfPtSd <- sd(dfPt$X.Pt.)

#durchschnittlicher Abstand zwischen den dfPts
dtdfPt <- t/length(dfPtnumeric)
#um die Zeit in das Koordinatensystem mit dfPt eintragen zu kˆnnen m¸ssen wir alle eintr‰ge * dtdfPt rechnen
i <- 1
zeitdfPt <- NULL
while (i <= length(dfPtnumeric)){
  h <- i*dtdfPt
  zeitdfPt <- rbind(zeitdfPt, data.frame(h))
  i <- i+1
}

#plot f¸r Schrittklassifizierung
plot(zeitdfPt$h, dfPtnumeric, type="p", col="black", xlab = "Zeit", ylab = "Schrittl‰nge aus Integration")

#Schrittklassifizierung
i <- 1
limit <- length(dfPtnumeric)
if((dfPtSd/dfPtmean)>(1/3)) {
  abline(h=dfPtmean+(dfPtSd/2), col="orange")
  abline(h=dfPtmean-(dfPtSd/2), col="green")
  abline(h=dfPtmean, col="blue")
  schrittegross <- c()
  schritteklein <- c()
  schrittenormal <- c()
  while(i<=limit) {
    if (dfPtnumeric[i]>(dfPtmean+(dfPtSd/2))){
      cat(i,dfPtnumeric[i],sep = '  ',"Groﬂer Schritt\n")
      schrittegross <- c(schrittegross, dfPtnumeric[i])
    }else if(dfPtnumeric[i]<(dfPtmean-(dfPtSd/2))){
      cat(i,dfPtnumeric[i],sep = '  ',"kleiner Schritt\n")
      schritteklein <- c(schritteklein, dfPtnumeric[i])
    }else{
      cat(i,dfPtnumeric[i],sep = '  ',"normaler Schritt\n")
      schrittenormal <- c(schrittenormal, dfPtnumeric[i])
    }
    i <- i+1
  }
  cat("\n")
  cat(length(schrittegross), "Groﬂe Schritte\n" , sep=' ')
  cat(length(schritteklein), "Kleine Schritte\n" , sep=' ')
  cat(length(schrittenormal), "Normale Schritte\n" , sep=' ')
}else {
  print("Schritte sind alle ‰hnlich groﬂ")
}


numSteps <- trunc(length(x.points)/2)
StreckeIntegration <- sum(dfPtnumeric)


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
  if((abs(dfPtnumeric[x]-dfPtnumeric[y])>0.1) & (abs(dfPtnumeric[x]-dfPtnumeric[j])>0.1) | (y == numSteps)) {
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

#Schritte insgesamt
cat(numSteps, sep = '', " Schritte insgesamt")

#Strecke insgesamt (mit Integration berechnet)
cat(round(StreckeIntegration, digits = 3),"m", sep = '', " Strecke errechnet mit Integration")
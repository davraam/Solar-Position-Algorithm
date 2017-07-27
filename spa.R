#' 
#' @title spa (Solar Position Algorithm)
#' @description The function implements the Solar Position Algorithm for Solar Radiation applications as 
#' described as a step by step procedure in Ibrahim R and Afshin A, (2008) technical report.
#' @details The Solar Position Algorithm calculates the solar zenith and azimuth angle with uncertainties
#' equal to $\pm 0.0003$ degrees in the period from the year -2000 to 6000. In this algorithm, the azimuth angle
#' is measured eastward from north and the observer's geographical longitude is considered negative west, or
#' positive east from Greenwich.
#' @param timestamp, a POSIXct object inticating the observer local time and date.
#' @param timezone, observer time zone, integer valid in the range -18 to 18 hours.
#' Negative values indicate zones west of Greenwich. 
#' @param longitude, the observer longitude in degrees valid in the range -180 to 180 degrees. 
#' Negative values indicate location west of Greenwich.
#' @param latitude, the observer latitude in degrees valid in the range -90 to 90 degrees. 
#' Negative values indicate location south of Equator.
#' @param elevation, the observer elevation in meters.
#' @param pressure, the annual average local pressure in millibars.
#' @param temperature, the annual average local temperature in Celcius degrees.
#' @param surface.slope, the Surface slope in degrees (measured from the horizontal plane), valid in the range
#'  -360 to 360 degrees.
#' @param surface.azimuth.rotation, the Surface azimuth rotation in degrees (measured from south to projection
#' of surface normal on horizontal plane), valid in the range -360 to 360 degrees.
#' @param DeltaT, in seconds. Is the difference between the Earth rotation time and Terrestrial Time.
#' It is derived from observation only and reported yearly in the Astronomical Almanac 
#' @return a list with (i) the topocentric zenith angle, theta (in degrees), (ii) the topocentric azimuth angle,
#' Phi (in degrees), (iii) the incidence angle for a surface oriented in any direction, I (in degrees), (iv) the
#' sunrise, R (in local time), and (v) the sunset, S (in local time).
#' @author Avraam, D. 
#' @export
#' @examples {
#'
#'   # Example 1
#'   timestamp <- as.POSIXct("2003-10-17 12:30:30")
#'   spa(timestamp, timezone=-7, longitude=-105.1786, latitude=39.742476, elevation=1830.14, pressure=820, temperature=11, surface.slope=30, surface.azimuth.rotation=-10, DeltaT=67)
#'
#' }
#'
spa <- function(timestamp=NULL, timezone=NULL, longitude=NULL, latitude=NULL, elevation=NULL, pressure=NULL, temperature=NULL, surface.slope=NULL, surface.azimuth.rotation=NULL, DeltaT=NULL){

# Constant parameters
DeltaUT1 <- 0
B <- 0
SunRadius <- 0.26667

# Extract year, month, day, hout, minute, seconds and timezone from the input timestamp arqument which is a POSIXct object 
year <- as.numeric(strftime(timestamp, format="%Y"))
month <- as.numeric(strftime(timestamp, format="%m"))
day <- as.numeric(strftime(timestamp, format="%d"))
hour <- as.numeric(strftime(timestamp, format="%H"))
minute <- as.numeric(strftime(timestamp, format="%M"))
sec <- as.numeric(strftime(timestamp, format="%S"))
#timezone <- as.numeric(strftime(timestamp, format="%z"))/100  # need change!

tz <- timezone

# Input parameters
Y <- year
M <- month
P <- pressure  # the annual average local pressure (in millibars)
T <- temperature  # the annual average local temperature (in degrees C)
omega <- surface.slope  # the slope of the surface measured from the horizontal plane
gamma <- surface.azimuth.rotation  # the surface azimuth rotation angle, measured from south to the projection of the surface

# Option to return all the values with 20 decimal places
options(digits=20)


########################################################################
# 
# Calculate the Julian and Julian Ephemeris Day, Century and Millennium:
#
########################################################################

# convert day to day in decimals (D)
D <- day + (hour - tz + (minute + (sec + DeltaUT1)/60)/60)/24

if(month < 3){
  M <- month + 12
  Y <- year - 1
}

# Calcuate the Julian Day (JD)
JD <-  trunc(365.25*(Y+4716)) + trunc(30.6001*(M+1)) + D + B - 1524.5

if (JD > 2299160){
  A <- trunc(Y/100)
  JD <- JD + (2 - A + trunc(A/4))
}

# Calculate the Julian Ephemeris Day (JDE)
JDE <- JD + DeltaT/86400

# Calculate the Julian Century (JC) and the Julian Ephemeris Century (JCE) for the 2000 standard epoch
JC <- (JD - 2451545)/36525
JCE <- (JDE - 2451545)/36525

# Calculate the Julian Ephemeris Millenium (JME) for the 2000 standard epoch
JME <- JCE/10


########################################################################
# 
# Calculate the Earth heliocentric longitude, latitude, and radius 
# vector (L, B, and R):
#
########################################################################

# Read the values of Earth Periodic Terms (EPT) from a table uploaded on Github 
url1 <- 'https://raw.github.com/davraam/Solar-Position-Algorithm/master/earth.periodic.terms.csv'
EPT <- read.csv(url1)

# Calculate the term L0 (in radians)
EPT_L0 <- EPT[which(EPT$Term=='L0'),]  # create a subset with L0 terms
L0 <- c()
for (i in 1:dim(EPT_L0)[1]){
  L0[i] <- EPT_L0[i,2] * cos(EPT_L0[i,3]+EPT_L0[i,4]*JME) 
}
L0 <- sum(L0)

# Calculate the term L1 (in radians)
EPT_L1 <- EPT[which(EPT$Term=='L1'),]  # create a subset with L1 terms
L1 <- c()
for (i in 1:dim(EPT_L1)[1]){
  L1[i] <- EPT_L1[i,2] * cos(EPT_L1[i,3]+EPT_L1[i,4]*JME) 
}
L1 <- sum(L1)

# Calculate the term L2 (in radians)
EPT_L2 <- EPT[which(EPT$Term=='L2'),]  # create a subset with L2 terms
L2 <- c()
for (i in 1:dim(EPT_L2)[1]){
  L2[i] <- EPT_L2[i,2] * cos(EPT_L2[i,3]+EPT_L2[i,4]*JME) 
}
L2 <- sum(L2)

# Calculate the term L3 (in radians)
EPT_L3 <- EPT[which(EPT$Term=='L3'),]  # create a subset with L3 terms
L3 <- c()
for (i in 1:dim(EPT_L3)[1]){
  L3[i] <- EPT_L3[i,2] * cos(EPT_L3[i,3]+EPT_L3[i,4]*JME) 
}
L3 <- sum(L3)

# Calculate the term L4 (in radians)
EPT_L4 <- EPT[which(EPT$Term=='L4'),]  # create a subset with L4 terms
L4 <- c()
for (i in 1:dim(EPT_L4)[1]){
  L4[i] <- EPT_L4[i,2] * cos(EPT_L4[i,3]+EPT_L4[i,4]*JME) 
}
L4 <- sum(L4)

# Calculate the term L5 (in radians)
EPT_L5 <- EPT[which(EPT$Term=='L5'),]  # create a subset with L5 terms
L5 <- c()
for (i in 1:dim(EPT_L5)[1]){
  L5[i] <- EPT_L5[i,2] * cos(EPT_L5[i,3]+EPT_L5[i,4]*JME) 
}
L5 <- sum(L5)

# Calculate the Earth heliocentric longitude, L (in radians)
L <- (L0+L1*JME+L2*JME^2+L3*JME^3+L4*JME^4+L5*JME^5)/(10^8)

# Calculate L in degrees
L <- L*180/pi

# Limit L to the range from 0 to 360 degrees
F <- L/360 - trunc(L/360)
if (L<0){
  L.limited <- 360-360*abs(F)
}else{
  L.limited <- 360*F
}

# Calculate the term B0 (in radians)
EPT_B0 <- EPT[which(EPT$Term=='B0'),]  # create a subset with B0 terms
B0 <- c()
for (i in 1:dim(EPT_B0)[1]){
  B0[i] <- EPT_B0[i,2] * cos(EPT_B0[i,3]+EPT_B0[i,4]*JME) 
}
B0 <- sum(B0)

# Calculate the term B1 (in radians)
EPT_B1 <- EPT[which(EPT$Term=='B1'),]  # create a subset with B1 terms
B1 <- c()
for (i in 1:dim(EPT_B1)[1]){
  B1[i] <- EPT_B1[i,2] * cos(EPT_B1[i,3]+EPT_B1[i,4]*JME) 
}
B1 <- sum(B1)

# Calculate the Earth heliocentric latitude, B (in degrees)
B <- (B0+B1*JME)/(10^8)

# Calculate B in degrees
B <- B*180/pi

# Calculate the term R0
EPT_R0 <- EPT[which(EPT$Term=='R0'),]  # create a subset with R0 terms
R0 <- c()
for (i in 1:dim(EPT_R0)[1]){
  R0[i] <- EPT_R0[i,2] * cos(EPT_R0[i,3]+EPT_R0[i,4]*JME) 
}
R0 <- sum(R0)

# Calculate the term R1
EPT_R1 <- EPT[which(EPT$Term=='R1'),]  # create a subset with R1 terms
R1 <- c()
for (i in 1:dim(EPT_R1)[1]){
  R1[i] <- EPT_R1[i,2] * cos(EPT_R1[i,3]+EPT_R1[i,4]*JME) 
}
R1 <- sum(R1)

# Calculate the term R2
EPT_R2 <- EPT[which(EPT$Term=='R2'),]  # create a subset with R2 terms
R2 <- c()
for (i in 1:dim(EPT_R2)[1]){
  R2[i] <- EPT_R2[i,2] * cos(EPT_R2[i,3]+EPT_R2[i,4]*JME) 
}
R2 <- sum(R2)

# Calculate the term R3
EPT_R3 <- EPT[which(EPT$Term=='R3'),]  # create a subset with R3 terms
R3 <- c()
for (i in 1:dim(EPT_R3)[1]){
  R3[i] <- EPT_R3[i,2] * cos(EPT_R3[i,3]+EPT_R3[i,4]*JME) 
}
R3 <- sum(R3)

# Calculate the term R4
EPT_R4 <- EPT[which(EPT$Term=='R4'),]  # create a subset with R4 terms
R4 <- c()
for (i in 1:dim(EPT_R4)[1]){
  R4[i] <- EPT_R4[i,2] * cos(EPT_R4[i,3]+EPT_R4[i,4]*JME) 
}
R4 <- sum(R4)

# Calculate the Earth radius vector, R (in Astronomical Units, AU)
R <- (R0+R1*JME+R2*JME^2+R3*JME^3+R4*JME^4)/(10^8)


########################################################################
# 
# Calculate the geocentric longitude and latitude (Theta and beta):
#
########################################################################

# Calculate the geocentric longitude, Theta (in degrees)
Theta <- L.limited + 180

# Limit Theta to the range from 0 to 360 degrees
if (Theta<0){
  Theta.limited <- 360 - 360*(abs(Theta/360 - trunc(Theta/360)))
}else{
  Theta.limited <- 360*(Theta/360 - trunc(Theta/360))
}

# Calculate the geocentric latitude, beta (in degrees)
beta <- -B


########################################################################
# 
# Calculate the nutation in longitude and obliquity (Deltapsi and 
# Deltaepsilon):
#
########################################################################

# Calculate the mean elongation of the moon from the sun, X0 (in degrees)
X0 <- 297.85036 + 445267.111480*JCE - 0.0019142*(JCE^2) + (JCE^3)/189474

# Calculate the mean anomaly of the sun (Earth), X1 (in degrees)
X1 <- 357.52772 + 35999.050340*JCE - 0.0001603*(JCE^2) - (JCE^3)/300000

# Calculate the mean anomaly of the moon, X2 (in degrees)
X2 <- 134.96298 + 477198.867398*JCE + 0.0086972*(JCE^2) + (JCE^3)/56250

# Calculate the moon's argument of latitude, X3 (in degrees)
X3 <- 93.27191 + 483202.017538*JCE - 0.0036825*(JCE^2) + (JCE^3)/327270

# Calculate the longitude of the ascending node of the moon's mean orbit on the ecliptic,
# measured from the mean equinox of the date, X4 (in degrees)
X4 <- 125.04452 - 1934.136261*JCE + 0.0020708*(JCE^2) + (JCE^3)/450000

# Read the values of Periodic Terms for the Nutation in Longitude and Obliquity (NLO) from a table uploaded on Github 
url2 <- 'https://raw.github.com/davraam/Solar-Position-Algorithm/master/nutation.longitude.obliquity.csv'
NLO <- read.csv(url2)

# Calculate the terms Deltapsi (in 0.0001 of arc seconds)
Deltapsi <- c()
for (i in 1:dim(NLO)[1]){
  Deltapsi[i] <- (NLO[i,6] + NLO[i,7]*JCE) * sin(pi*(X0*NLO[i,1]+X1*NLO[i,2]+X2*NLO[i,3]+X3*NLO[i,4]+X4*NLO[i,5])/180)
}

# Calculate the nutation in longitude, Deltapsi (in degrees)
Deltapsi <- sum(Deltapsi)/36000000

# Calculate the terms Deltaepsilon (in 0.0001 of arc seconds)
Deltaepsilon <- c()
for (i in 1:dim(NLO)[1]){
  Deltaepsilon[i] <- (NLO[i,8] + NLO[i,9]*JCE) * cos(pi*(X0*NLO[i,1]+X1*NLO[i,2]+X2*NLO[i,3]+X3*NLO[i,4]+X4*NLO[i,5])/180)
}

# Calculate the nutation in obliquity, Deltaepsilon (in degrees)
Deltaepsilon <- sum(Deltaepsilon)/36000000


########################################################################
# 
# Calculate the true obliquity of the ecliptic, epsilon (in degrees):
#
########################################################################

# Calculate the true obliquity of the ecliptic, epsilon0 (in arc seconds)
U <- JME/10
epsilon0 <- 84381.448 - 4680.93*U - 1.55*U^2 + 1999.25*U^3 - 51.38*U^4 -249.67*U^5 - 39.05*U^6 + 7.12*U^7 + 27.87*U^8 + 5.79*U^9 + 2.45*U^10

# Calculate the true obliquity of the ecliptic, epsilon (in degrees) 
epsilon <- epsilon0/3600 + Deltaepsilon


########################################################################
# 
# Calculate the aberration correction, Deltatau (in degrees):
#
########################################################################

Deltatau <- -20.4898/(3600*R)


########################################################################
# 
# Calculate the apparent sun longitude, lambda (in degrees):
#
########################################################################

lambda <- Theta.limited + Deltapsi + Deltatau


########################################################################
# 
# Calculate the apparent sidereal time at Greenwich at any given time, v
# (in degrees):
#
########################################################################

# Calculate the mean sidereal time at Greenwich, v0 (in degrees)
v0 <- 280.46061837 + 360.98564736629*(JD-2451545) + 0.000387933*JC^2 - (JC^3)/38710000

# Limit v0 to the range from 0 to 360 degrees
if (v0 < 0){
  v0.limited <- 360 - 360*(abs(v0/360 - trunc(v0/360)))
}else{
  v0.limited <- 360*(v0/360 - trunc(v0/360))
}

# Calculate the apparent sidereal time in Greenwich, v (in degrees)
v <- v0.limited + Deltapsi*cos(epsilon*pi/180)


########################################################################
# 
# Calculate the geocentric sun right ascension, alpha (in degrees):
#
########################################################################

# Calculate the sun right ascension, alpha (in radians)
alpha <- atan2((sin(lambda*pi/180)*cos(epsilon*pi/180)-tan(beta*pi/180)*sin(epsilon*pi/180)),cos(lambda*pi/180))

# Calculate alpha in degrees and limit it to the range from 0 to 360 degrees
alpha <- (alpha*180)/pi
if (alpha<0){
  alpha.limited <- 360 - 360*(abs(alpha/360 - trunc(alpha/360)))
}else{
  alpha.limited <- 360*(alpha/360 - trunc(alpha/360))
}


########################################################################
# 
# Calculate the geocentric sun declination, delta (in degrees):
#
########################################################################

delta <- asin(sin(beta*pi/180)*cos(epsilon*pi/180)+cos(beta*pi/180)*sin(epsilon*pi/180)*sin(lambda*pi/180))
delta <- (delta*180)/pi


########################################################################
# 
# Calculate the observer local hout angle, H (in degrees):
#
########################################################################

sigma <- longitude
H <- v + sigma - alpha

# Limit H to the range from 0 to 360 degrees
if (H < 0){
  H.limited <- 360 - 360*(abs(H/360 - trunc(H/360)))
}else{
  H.limited <- 360*(H/360 - trunc(H/360))
}


########################################################################
# 
# Calculate the topocentric sun right ascension, alphaPrime (in degrees):
#
########################################################################

# Calculate the equatorial horizontal parallax of the sun, xi (in degrees)
xi <- 8.794/(3600*R)

# Calculate the term u (in radians)
phi <- latitude
u <- atan(0.99664719*tan(phi*pi/180))

# Calculate the term x
E <- elevation
x <- cos(u) + (E/6378140)*cos(phi*pi/180)

# Calculate the term y
y <- 0.99664719*sin(u) + (E/6378140)*sin(phi*pi/180)

# Calculate the parallax in sun right ascension, Deltaalpha (in degrees)
Deltaalpha <- atan2((-x*sin(xi*pi/180)*sin(H*pi/180)),(cos(delta*pi/180)-x*sin(xi*pi/180)*cos(H*pi/180)))
Deltaalpha <- Deltaalpha*180/pi

# Calculate the topocentric sun right ascension, alphaPrime (in degrees)
alphaPrime <- alpha.limited + Deltaalpha

# Calculate the topocentric sun declination, deltaPrime (in degrees)
deltaPrime <- atan2(((sin(delta*pi/180)-y*sin(xi*pi/180))*cos(Deltaalpha*pi/180)),(cos(delta*pi/180)-x*sin(xi*pi/180)*cos(H*pi/180)))
deltaPrime <- deltaPrime*180/pi


########################################################################
# 
# Calculate the topocentric local hour angle, HPrime (in degrees):
#
########################################################################
HPrime <- H.limited - Deltaalpha


########################################################################
# 
# Calculate the topocentric zenith angle, theta (in degrees):
#
########################################################################

# Calcuate the topocentric elevation angle without atmospheric refraction correction, e0 (in degrees)
e0 <- asin(sin(phi*pi/180)*sin(deltaPrime*pi/180)+cos(phi*pi/180)*cos(deltaPrime*pi/180)*cos(HPrime*pi/180))
e0 <- e0*180/pi
 
# Calculate the atmospheric refraction correction, Deltae (in degrees)
Deltae <- (P/1010)*(283/(273+T))*(1.02/(60*tan((e0+10.3/(e0+5.11))*pi/180)))

# Calcuate the topocentric elevation angle, e (in degrees)
e <- e0 + Deltae

# Calculate the topocentric zenith angle, theta (in degrees)
theta <- 90 - e


########################################################################
# 
# Calculate the topocentric azimuth angle, Phi (in degrees):
#
########################################################################

# Calculate the topocentric astronomers azimuth angle, Gamma (in degrees)
Gamma <- atan2(sin(HPrime*pi/180), (cos(HPrime*pi/180)*sin(phi*pi/180)-tan(deltaPrime*pi/180)*cos(phi*pi/180)))
Gamma <- Gamma*180/pi

# Limit Gamma to the range from 0 to 360 degrees
if (Gamma < 0){
  Gamma.limited <- 360 - 360*(abs(Gamma/360 - trunc(Gamma/360)))
}else{
  Gamma.limited <- 360*(Gamma/360 - trunc(Gamma/360))
}

# Calculate the topocentric azimuth angle, Phi for navigators and solar radiation users (in degrees)
Phi <- Gamma.limited + 180

# Limit Phi to the range from 0 to 360 degrees
if (Phi < 0){
  Phi.limited <- 360 - 360*(abs(Phi/360 - trunc(Phi/360)))
}else{
  Phi.limited <- 360*(Phi/360 - trunc(Phi/360))
}


########################################################################
# 
# Calculate the incidence angle for a surface oriented in any direction, 
# I (in degrees):
#
########################################################################

I <- acos(cos(theta*pi/180)*cos(omega*pi/180)+sin(omega*pi/180)*sin(theta*pi/180)*cos((Gamma.limited-gamma)*pi/180))
I <- I*180/pi


########################################################################
# 
# Equation of Time
#
########################################################################

# Calculate the sun's mean longitude, M (in degrees)
M <- 280.4664567 + 360007.6982779*JME + 0.03032028*JME^2 + (JME^3)/49931 - (JME^4)/15300 - (JME^5)/2000000

# Limit M to the range from 0 to 360 degrees
if (M < 0){
  M.limited <- 360 - 360*(abs(M/360 - trunc(M/360)))
}else{
  M.limited <- 360*(M/360 - trunc(M/360))
}

# Calculate the equation of time, E (in degrees) (E is the difference between solar apparent and mean time)
E <- M.limited - 0.0057183 - alpha + Deltapsi*cos(epsilon*pi/180)

# Multiply E by 4 to change its units from degrees to minutes of time
E <- E*4 

# Limit E if its absolute value is greater than 20
if (E < -20){ E <- E + 1440}
if (E > 20){ E <- E - 1440}


########################################################################
# 
# Sunrise, Sun Transit, and Sunset 
#
########################################################################


###############################################################################
#
# Calculate the apparent sidereal time at Greenwich at 0 UT, v.UT0 (in degrees)
# The code below repeats the steps 3.1.1-3.1.4, 3.4.1-3.5.2 and 3.8.1-3.8.3 
# from Ibrahim and Afshin 2008 Report. All the variables below have the suffix 
# .UT0
#
###############################################################################

# At 0 UT
DeltaT.UT0 <- 0 
D.UT0 <- day 
if(month < 3){
  M.UT0 <- month + 12
  Y.UT0 <- year - 1
}else{
  M.UT0 <- month
  Y.UT0 <- year
}

# Calcuate the Julian Day (JD.UT0)
JD.UT0 <-  trunc(365.25*(Y.UT0+4716)) + trunc(30.6001*(M.UT0+1)) + D.UT0 + B - 1524.5

if (JD.UT0 > 2299160){
  A.UT0 <- trunc(Y.UT0/100)
  JD.UT0 <- JD.UT0 + (2 - A.UT0 + trunc(A.UT0/4))
}

# Calculate the Julian Ephemeris Day (JDE.UT0)
JDE.UT0 <- JD.UT0 + DeltaT.UT0/86400

# Calculate the Julian Century (JC.UT0) and the Julian Ephemeris Century (JCE.UT0) for the 2000 standard epoch
JC.UT0 <- (JD.UT0 - 2451545)/36525
JCE.UT0 <- (JDE.UT0 - 2451545)/36525

# Calculate the Julian Ephemeris Millenium (JME.UT0) for the 2000 standard epoch
JME.UT0 <- JCE.UT0/10

# Calculate the mean elongation of the moon from the sun, X0.UT0 (in degrees)
X0.UT0 <- 297.85036 + 445267.111480*JCE.UT0 - 0.0019142*(JCE.UT0^2) + (JCE.UT0^3)/189474

# Calculate the mean anomaly of the sun (Earth), X1.UT0 (in degrees)
X1.UT0 <- 357.52772 + 35999.050340*JCE.UT0 - 0.0001603*(JCE.UT0^2) - (JCE.UT0^3)/300000

# Calculate the mean anomaly of the moon, X2.UT0 (in degrees)
X2.UT0 <- 134.96298 + 477198.867398*JCE.UT0 + 0.0086972*(JCE.UT0^2) + (JCE.UT0^3)/56250

# Calculate the moon's argument of latitude, X3.UT0 (in degrees)
X3.UT0 <- 93.27191 + 483202.017538*JCE.UT0 - 0.0036825*(JCE.UT0^2) + (JCE.UT0^3)/327270

# Calculate the longitude of the ascending node of the moon's mean orbit on the ecliptic,
# measured from the mean equinox of the date, X4.UT0 (in degrees)
X4.UT0 <- 125.04452 - 1934.136261*JCE.UT0 + 0.0020708*(JCE.UT0^2) + (JCE.UT0^3)/450000

# Calculate the terms Deltapsi.UT0 (in 0.0001 of arc seconds)
Deltapsi.UT0 <- c()
for (i in 1:dim(NLO)[1]){
  Deltapsi.UT0[i] <- (NLO[i,6] + NLO[i,7]*JCE.UT0) * sin(pi*(X0.UT0*NLO[i,1]+X1.UT0*NLO[i,2]+X2.UT0*NLO[i,3]+X3.UT0*NLO[i,4]+X4.UT0*NLO[i,5])/180)
}

# Calculate the nutation in longitude, Deltapsi.UT0 (in degrees)
Deltapsi.UT0 <- sum(Deltapsi.UT0)/36000000

# Calculate the terms Deltaepsilon.UT0 (in 0.0001 of arc seconds)
Deltaepsilon.UT0 <- c()
for (i in 1:dim(NLO)[1]){
  Deltaepsilon.UT0[i] <- (NLO[i,8] + NLO[i,9]*JCE.UT0) * cos(pi*(X0.UT0*NLO[i,1]+X1.UT0*NLO[i,2]+X2.UT0*NLO[i,3]+X3.UT0*NLO[i,4]+X4.UT0*NLO[i,5])/180)
}

# Calculate the nutation in obliquity, Deltaepsilon.UT0 (in degrees)
Deltaepsilon.UT0 <- sum(Deltaepsilon.UT0)/36000000

# Calculate the true obliquity of the ecliptic, epsilon0.UT0 (in arc seconds)
U.UT0 <- JME.UT0/10
epsilon0.UT0 <- 84381.448 - 4680.93*U.UT0 - 1.55*U.UT0^2 + 1999.25*U.UT0^3 - 51.38*U.UT0^4 -249.67*U.UT0^5 - 39.05*U.UT0^6 + 7.12*U.UT0^7 + 27.87*U.UT0^8 + 5.79*U.UT0^9 + 2.45*U.UT0^10

# Calculate the true obliquity of the ecliptic, epsilon (in degrees) 
epsilon.UT0 <- epsilon0.UT0/3600 + Deltaepsilon.UT0

# Calculate the mean sidereal time at Greenwich, v0.UT0 (in degrees)
v0.UT0 <- 280.46061837 + 360.98564736629*(JD.UT0-2451545) + 0.000387933*JC.UT0^2 - (JC.UT0^3)/38710000

# Limit v0.UT0 to the range from 0 to 360 degrees
if (v0.UT0 < 0){
  v0.UT0.limited <- 360 - 360*(abs(v0.UT0/360 - trunc(v0.UT0/360)))
}else{
  v0.UT0.limited <- 360*(v0.UT0/360 - trunc(v0.UT0/360))
}

# Calculate the apparent sidereal time in Greenwich, v.UT0 (in degrees)
v.UT0 <- v0.UT0.limited + Deltapsi.UT0*cos(epsilon.UT0*pi/180)

###############################################################################




###############################################################################
#
# Calculate the geocentric right ascension and declination at 0 TT for the day
# before the day of interest, the day of interest and the day after the day of
# interest. The code below repeats equations 3.1.1-3.7 and 3.9.1-3.10 from 
# Ibrahim and Afshin 2008 Report. All the variables below have the suffix 
# .TT0
#
###############################################################################

DeltaT.TT0 <- DeltaT

# Find the date before and after the day of interest
date <- as.Date(paste0(year, '-', month, '-', day))  # input date
date.before <- date-1  # the day before the day of interest
date.after <- date+1  # the day after the day of interest

D.bef.TT0 <- as.numeric(strftime(date.before, format="%d"))
M.bef.TT0 <- as.numeric(strftime(date.before, format="%m"))
Y.bef.TT0 <- as.numeric(strftime(date.before, format="%Y"))

if(M.bef.TT0 < 3){
  M.bef.TT0 <- M.bef.TT0 + 12
  Y.bef.TT0 <- Y.bef.TT0 - 1
}else{
  M.bef.TT0 <- M.bef.TT0
  Y.bef.TT0 <- Y.bef.TT0
}

D.aft.TT0 <- as.numeric(strftime(date.after, format="%d"))
M.aft.TT0 <- as.numeric(strftime(date.after, format="%m"))
Y.aft.TT0 <- as.numeric(strftime(date.after, format="%Y"))

if(M.aft.TT0 < 3){
  M.aft.TT0 <- M.aft.TT0 + 12
  Y.aft.TT0 <- Y.aft.TT0 - 1
}else{
  M.aft.TT0 <- M.aft.TT0
  Y.aft.TT0 <- Y.aft.TT0
}

D.TT0 <- as.numeric(strftime(date, format="%d"))
M.TT0 <- as.numeric(strftime(date, format="%m"))
Y.TT0 <- as.numeric(strftime(date, format="%Y"))

if(M.TT0 < 3){
  M.TT0 <- M.TT0 + 12
  Y.TT0 <- Y.TT0 - 1
}else{
  M.TT0 <- M.TT0
  Y.TT0 <- Y.TT0
}


# Calcuate the Julian Days (JD.bef.TT0, JD.TT0, JD.aft.TT0)
JD.bef.TT0 <-  trunc(365.25*(Y.bef.TT0+4716)) + trunc(30.6001*(M.bef.TT0+1)) + D.bef.TT0 + B - 1524.5
JD.TT0 <-  trunc(365.25*(Y.TT0+4716)) + trunc(30.6001*(M.TT0+1)) + D.TT0 + B - 1524.5
JD.aft.TT0 <-  trunc(365.25*(Y.aft.TT0+4716)) + trunc(30.6001*(M.aft.TT0+1)) + D.aft.TT0 + B - 1524.5

if (JD.bef.TT0 > 2299160){
  A.bef.TT0 <- trunc(Y.bef.TT0/100)
  JD.bef.TT0 <- JD.bef.TT0 + (2 - A.bef.TT0 + trunc(A.bef.TT0/4))
}

if (JD.TT0 > 2299160){
  A.TT0 <- trunc(Y.TT0/100)
  JD.TT0 <- JD.TT0 + (2 - A.TT0 + trunc(A.TT0/4))
}

if (JD.aft.TT0 > 2299160){
  A.aft.TT0 <- trunc(Y.aft.TT0/100)
  JD.aft.TT0 <- JD.aft.TT0 + (2 - A.aft.TT0 + trunc(A.aft.TT0/4))
}

# Calculate the Julian Ephemeris Days (JDE.bef.TT0, JDE.TT0, JDE.aft.TT0)
JDE.bef.TT0 <- JD.bef.TT0 + DeltaT.TT0/86400
JDE.TT0 <- JD.TT0 + DeltaT.TT0/86400
JDE.aft.TT0 <- JD.aft.TT0 + DeltaT.TT0/86400

# Calculate the Julian Century (JC.bef.TT0, JC.TT0, JC.aft.TT0) and the Julian Ephemeris Century (JCE.bef.TT0, JCE.TT0, JCE.aft.TT0) for the 2000 standard epoch
JC.bef.TT0 <- (JD.bef.TT0 - 2451545)/36525
JCE.bef.TT0 <- (JDE.bef.TT0 - 2451545)/36525
JC.TT0 <- (JD.TT0 - 2451545)/36525
JCE.TT0 <- (JDE.TT0 - 2451545)/36525
JC.aft.TT0 <- (JD.aft.TT0 - 2451545)/36525
JCE.aft.TT0 <- (JDE.aft.TT0 - 2451545)/36525

# Calculate the Julian Ephemeris Millenium (JME.TT0) for the 2000 standard epoch
JME.bef.TT0 <- JCE.bef.TT0/10
JME.TT0 <- JCE.TT0/10
JME.aft.TT0 <- JCE.aft.TT0/10

# Calculate the terms L0.bef.TT0, L0.TT0, L0.aft.TT0 (in radians)
L0.bef.TT0 <- c()
L0.TT0 <- c()
L0.aft.TT0 <- c()
for (i in 1:dim(EPT_L0)[1]){
  L0.bef.TT0[i] <- EPT_L0[i,2] * cos(EPT_L0[i,3]+EPT_L0[i,4]*JME.bef.TT0)
  L0.TT0[i] <- EPT_L0[i,2] * cos(EPT_L0[i,3]+EPT_L0[i,4]*JME.TT0) 
  L0.aft.TT0[i] <- EPT_L0[i,2] * cos(EPT_L0[i,3]+EPT_L0[i,4]*JME.aft.TT0)  
}
L0.bef.TT0 <- sum(L0.bef.TT0)
L0.TT0 <- sum(L0.TT0)
L0.aft.TT0 <- sum(L0.aft.TT0)

# Calculate the terms L1.bef.TT0, L1.TT0, L1.aft.TT0 (in radians)
L1.bef.TT0 <- c()
L1.TT0 <- c()
L1.aft.TT0 <- c()
for (i in 1:dim(EPT_L1)[1]){
  L1.bef.TT0[i] <- EPT_L1[i,2] * cos(EPT_L1[i,3]+EPT_L1[i,4]*JME.bef.TT0)
  L1.TT0[i] <- EPT_L1[i,2] * cos(EPT_L1[i,3]+EPT_L1[i,4]*JME.TT0) 
  L1.aft.TT0[i] <- EPT_L1[i,2] * cos(EPT_L1[i,3]+EPT_L1[i,4]*JME.aft.TT0)  
}
L1.bef.TT0 <- sum(L1.bef.TT0)
L1.TT0 <- sum(L1.TT0)
L1.aft.TT0 <- sum(L1.aft.TT0)

# Calculate the terms L2.bef.TT0, L2.TT0, L2.aft.TT0 (in radians)
L2.bef.TT0 <- c()
L2.TT0 <- c()
L2.aft.TT0 <- c()
for (i in 1:dim(EPT_L2)[1]){
  L2.bef.TT0[i] <- EPT_L2[i,2] * cos(EPT_L2[i,3]+EPT_L2[i,4]*JME.bef.TT0)
  L2.TT0[i] <- EPT_L2[i,2] * cos(EPT_L2[i,3]+EPT_L2[i,4]*JME.TT0) 
  L2.aft.TT0[i] <- EPT_L2[i,2] * cos(EPT_L2[i,3]+EPT_L2[i,4]*JME.aft.TT0)  
}
L2.bef.TT0 <- sum(L2.bef.TT0)
L2.TT0 <- sum(L2.TT0)
L2.aft.TT0 <- sum(L2.aft.TT0)

# Calculate the terms L3.bef.TT0, L3.TT0, L3.aft.TT0 (in radians)
L3.bef.TT0 <- c()
L3.TT0 <- c()
L3.aft.TT0 <- c()
for (i in 1:dim(EPT_L3)[1]){
  L3.bef.TT0[i] <- EPT_L3[i,2] * cos(EPT_L3[i,3]+EPT_L3[i,4]*JME.bef.TT0)
  L3.TT0[i] <- EPT_L3[i,2] * cos(EPT_L3[i,3]+EPT_L3[i,4]*JME.TT0) 
  L3.aft.TT0[i] <- EPT_L3[i,2] * cos(EPT_L3[i,3]+EPT_L3[i,4]*JME.aft.TT0)   
}
L3.bef.TT0 <- sum(L3.bef.TT0)
L3.TT0 <- sum(L3.TT0)
L3.aft.TT0 <- sum(L3.aft.TT0)

# Calculate the terms L4.bef.TT0, L4.TT0, L4.aft.TT0 (in radians)
L4.bef.TT0 <- c()
L4.TT0 <- c()
L4.aft.TT0 <- c()
for (i in 1:dim(EPT_L4)[1]){
  L4.bef.TT0[i] <- EPT_L4[i,2] * cos(EPT_L4[i,3]+EPT_L4[i,4]*JME.bef.TT0)
  L4.TT0[i] <- EPT_L4[i,2] * cos(EPT_L4[i,3]+EPT_L4[i,4]*JME.TT0) 
  L4.aft.TT0[i] <- EPT_L4[i,2] * cos(EPT_L4[i,3]+EPT_L4[i,4]*JME.aft.TT0)   
}
L4.bef.TT0 <- sum(L4.bef.TT0)
L4.TT0 <- sum(L4.TT0)
L4.aft.TT0 <- sum(L4.aft.TT0)

# Calculate the terms L5.bef.TT0, L5.TT0, L5.aft.TT0 (in radians)
L5.bef.TT0 <- c()
L5.TT0 <- c()
L5.aft.TT0 <- c()
for (i in 1:dim(EPT_L5)[1]){
  L5.bef.TT0[i] <- EPT_L5[i,2] * cos(EPT_L5[i,3]+EPT_L5[i,4]*JME.bef.TT0)
  L5.TT0[i] <- EPT_L5[i,2] * cos(EPT_L5[i,3]+EPT_L5[i,4]*JME.TT0) 
  L5.aft.TT0[i] <- EPT_L5[i,2] * cos(EPT_L5[i,3]+EPT_L5[i,4]*JME.aft.TT0)   
}
L5.bef.TT0 <- sum(L5.bef.TT0)
L5.TT0 <- sum(L5.TT0)
L5.aft.TT0 <- sum(L5.aft.TT0)

# Calculate the Earth heliocentric longitude, L.bef.TT0, L.TT0, L.aft.TT0 (in radians)
L.bef.TT0 <- (L0.bef.TT0+L1.bef.TT0*JME.bef.TT0+L2.bef.TT0*JME.bef.TT0^2+L3.bef.TT0*JME.bef.TT0^3+L4.bef.TT0*JME.bef.TT0^4+L5.bef.TT0*JME.bef.TT0^5)/(10^8)
L.TT0 <- (L0.TT0+L1.TT0*JME.TT0+L2.TT0*JME.TT0^2+L3.TT0*JME.TT0^3+L4.TT0*JME.TT0^4+L5.TT0*JME.TT0^5)/(10^8)
L.aft.TT0 <- (L0.aft.TT0+L1.aft.TT0*JME.aft.TT0+L2.aft.TT0*JME.aft.TT0^2+L3.aft.TT0*JME.aft.TT0^3+L4.aft.TT0*JME.aft.TT0^4+L5.aft.TT0*JME.aft.TT0^5)/(10^8)

# Calculate L.bef.TT0, L.TT0, L.aft.TT0 in degrees
L.bef.TT0 <- L.bef.TT0*180/pi
L.TT0 <- L.TT0*180/pi
L.aft.TT0 <- L.aft.TT0*180/pi

# Limit L.bef.TT0, L.TT0, L.aft.TT0 to the range from 0 to 360 degrees
F.bef.TT0 <- L.bef.TT0/360 - trunc(L.bef.TT0/360)
F.TT0 <- L.TT0/360 - trunc(L.TT0/360)
F.aft.TT0 <- L.aft.TT0/360 - trunc(L.aft.TT0/360)
if (L.bef.TT0 < 0){
  L.limited.bef.TT0 <- 360-360*abs(F.bef.TT0)
}else{
  L.limited.bef.TT0 <- 360*F.bef.TT0
}
if (L.TT0 < 0){
  L.limited.TT0 <- 360-360*abs(F.TT0)
}else{
  L.limited.TT0 <- 360*F.TT0
}
if (L.aft.TT0 < 0){
  L.limited.aft.TT0 <- 360-360*abs(F.aft.TT0)
}else{
  L.limited.aft.TT0 <- 360*F.aft.TT0
}


# Calculate the terms B0.bef.TT0, B0.TT0, B0.aft.TT0 (in radians)
B0.bef.TT0 <- c()
B0.TT0 <- c()
B0.aft.TT0 <- c()
for (i in 1:dim(EPT_B0)[1]){
  B0.bef.TT0[i] <- EPT_B0[i,2] * cos(EPT_B0[i,3]+EPT_B0[i,4]*JME.bef.TT0)
  B0.TT0[i] <- EPT_B0[i,2] * cos(EPT_B0[i,3]+EPT_B0[i,4]*JME.TT0) 
  B0.aft.TT0[i] <- EPT_B0[i,2] * cos(EPT_B0[i,3]+EPT_B0[i,4]*JME.aft.TT0) 
}
B0.bef.TT0 <- sum(B0.bef.TT0)
B0.TT0 <- sum(B0.TT0)
B0.aft.TT0 <- sum(B0.aft.TT0)

# Calculate the terms B1.bef.TT0, B1.TT0, B1.aft.TT0 (in radians)
B1.bef.TT0 <- c()
B1.TT0 <- c()
B1.aft.TT0 <- c()
for (i in 1:dim(EPT_B1)[1]){
  B1.bef.TT0[i] <- EPT_B1[i,2] * cos(EPT_B1[i,3]+EPT_B1[i,4]*JME.bef.TT0)
  B1.TT0[i] <- EPT_B1[i,2] * cos(EPT_B1[i,3]+EPT_B1[i,4]*JME.TT0) 
  B1.aft.TT0[i] <- EPT_B1[i,2] * cos(EPT_B1[i,3]+EPT_B1[i,4]*JME.aft.TT0) 
}
B1.bef.TT0 <- sum(B1.bef.TT0)
B1.TT0 <- sum(B1.TT0)
B1.aft.TT0 <- sum(B1.aft.TT0)

# Calculate the Earth heliocentric latitude, B.bef.TT0, B.TT0, B.aft.TT0 (in degrees)
B.bef.TT0 <- (B0.bef.TT0+B1.bef.TT0*JME.bef.TT0)/(10^8)
B.TT0 <- (B0.TT0+B1.TT0*JME.TT0)/(10^8)
B.aft.TT0 <- (B0.aft.TT0+B1.aft.TT0*JME.aft.TT0)/(10^8)

# Calculate B.bef.TT0, B.TT0, B.aft.TT0 in degrees
B.bef.TT0 <- B.bef.TT0*180/pi
B.TT0 <- B.TT0*180/pi
B.aft.TT0 <- B.aft.TT0*180/pi

# Calculate the terms R0.bef.TT0, R0.TT0, R0.aft.TT0
R0.bef.TT0 <- c()
R0.TT0 <- c()
R0.aft.TT0 <- c()
for (i in 1:dim(EPT_R0)[1]){
  R0.bef.TT0[i] <- EPT_R0[i,2] * cos(EPT_R0[i,3]+EPT_R0[i,4]*JME.bef.TT0)
  R0.TT0[i] <- EPT_R0[i,2] * cos(EPT_R0[i,3]+EPT_R0[i,4]*JME.TT0) 
  R0.aft.TT0[i] <- EPT_R0[i,2] * cos(EPT_R0[i,3]+EPT_R0[i,4]*JME.aft.TT0)   
}
R0.bef.TT0 <- sum(R0.bef.TT0)
R0.TT0 <- sum(R0.TT0)
R0.aft.TT0 <- sum(R0.aft.TT0)

# Calculate the terms R1.bef.TT0, R1.TT0, R1.aft.TT0
R1.bef.TT0 <- c()
R1.TT0 <- c()
R1.aft.TT0 <- c()
for (i in 1:dim(EPT_R1)[1]){
  R1.bef.TT0[i] <- EPT_R1[i,2] * cos(EPT_R1[i,3]+EPT_R1[i,4]*JME.bef.TT0)
  R1.TT0[i] <- EPT_R1[i,2] * cos(EPT_R1[i,3]+EPT_R1[i,4]*JME.TT0) 
  R1.aft.TT0[i] <- EPT_R1[i,2] * cos(EPT_R1[i,3]+EPT_R1[i,4]*JME.aft.TT0)   
}
R1.bef.TT0 <- sum(R1.bef.TT0)
R1.TT0 <- sum(R1.TT0)
R1.aft.TT0 <- sum(R1.aft.TT0)

# Calculate the terms R2.bef.TT0, R2.TT0, R2.aft.TT0
R2.bef.TT0 <- c()
R2.TT0 <- c()
R2.aft.TT0 <- c()
for (i in 1:dim(EPT_R2)[1]){
  R2.bef.TT0[i] <- EPT_R2[i,2] * cos(EPT_R2[i,3]+EPT_R2[i,4]*JME.bef.TT0)
  R2.TT0[i] <- EPT_R2[i,2] * cos(EPT_R2[i,3]+EPT_R2[i,4]*JME.TT0) 
  R2.aft.TT0[i] <- EPT_R2[i,2] * cos(EPT_R2[i,3]+EPT_R2[i,4]*JME.aft.TT0)   
}
R2.bef.TT0 <- sum(R2.bef.TT0)
R2.TT0 <- sum(R2.TT0)
R2.aft.TT0 <- sum(R2.aft.TT0)

# Calculate the terms R3.bef.TT0, R3.TT0, R3.aft.TT0
R3.bef.TT0 <- c()
R3.TT0 <- c()
R3.aft.TT0 <- c()
for (i in 1:dim(EPT_R3)[1]){
  R3.bef.TT0[i] <- EPT_R3[i,2] * cos(EPT_R3[i,3]+EPT_R3[i,4]*JME.bef.TT0)
  R3.TT0[i] <- EPT_R3[i,2] * cos(EPT_R3[i,3]+EPT_R3[i,4]*JME.TT0) 
  R3.aft.TT0[i] <- EPT_R3[i,2] * cos(EPT_R3[i,3]+EPT_R3[i,4]*JME.aft.TT0)   
}
R3.bef.TT0 <- sum(R3.bef.TT0)
R3.TT0 <- sum(R3.TT0)
R3.aft.TT0 <- sum(R3.aft.TT0)

# Calculate the terms R4.bef.TT0, R4.TT0, R4.aft.TT0
R4.bef.TT0 <- c()
R4.TT0 <- c()
R4.aft.TT0 <- c()
for (i in 1:dim(EPT_R4)[1]){
  R4.bef.TT0[i] <- EPT_R4[i,2] * cos(EPT_R4[i,3]+EPT_R4[i,4]*JME.bef.TT0)
  R4.TT0[i] <- EPT_R4[i,2] * cos(EPT_R4[i,3]+EPT_R4[i,4]*JME.TT0) 
  R4.aft.TT0[i] <- EPT_R4[i,2] * cos(EPT_R4[i,3]+EPT_R4[i,4]*JME.aft.TT0)   
}
R4.bef.TT0 <- sum(R4.bef.TT0)
R4.TT0 <- sum(R4.TT0)
R4.aft.TT0 <- sum(R4.aft.TT0)

# Calculate the Earth radius vector, R.bef.TT0, R.TT0, R.aft.TT0 (in Astronomical Units, AU)
R.bef.TT0 <- (R0.bef.TT0+R1.bef.TT0*JME.bef.TT0+R2.bef.TT0*JME.bef.TT0^2+R3.bef.TT0*JME.bef.TT0^3+R4.bef.TT0*JME.bef.TT0^4)/(10^8)
R.TT0 <- (R0.TT0+R1.TT0*JME.TT0+R2.TT0*JME.TT0^2+R3.TT0*JME.TT0^3+R4.TT0*JME.TT0^4)/(10^8)
R.aft.TT0 <- (R0.aft.TT0+R1.aft.TT0*JME.aft.TT0+R2.aft.TT0*JME.aft.TT0^2+R3.aft.TT0*JME.aft.TT0^3+R4.aft.TT0*JME.aft.TT0^4)/(10^8)

# Calculate the geocentric longitude, Theta.bef.TT0, Theta.TT0, Theta.aft.TT0 (in degrees)
Theta.bef.TT0 <- L.limited.bef.TT0 + 180
Theta.TT0 <- L.limited.TT0 + 180
Theta.aft.TT0 <- L.limited.aft.TT0 + 180

# Limit Theta.bef.TT0, Theta.TT0, Theta.aft.TT0 to the range from 0 to 360 degrees
if (Theta.bef.TT0 < 0){
  Theta.limited.bef.TT0 <- 360 - 360*(abs(Theta.bef.TT0/360 - trunc(Theta.bef.TT0/360)))
}else{
  Theta.limited.bef.TT0 <- 360*(Theta.bef.TT0/360 - trunc(Theta.bef.TT0/360))
}
if (Theta.TT0 < 0){
  Theta.limited.TT0 <- 360 - 360*(abs(Theta.TT0/360 - trunc(Theta.TT0/360)))
}else{
  Theta.limited.TT0 <- 360*(Theta.TT0/360 - trunc(Theta.TT0/360))
}
if (Theta.aft.TT0 < 0){
  Theta.limited.aft.TT0 <- 360 - 360*(abs(Theta.aft.TT0/360 - trunc(Theta.aft.TT0/360)))
}else{
  Theta.limited.aft.TT0 <- 360*(Theta.aft.TT0/360 - trunc(Theta.aft.TT0/360))
}

# Calculate the geocentric latitude, beta.bef.TT0, beta.TT0, beta.aft.TT0 (in degrees)
beta.bef.TT0 <- -B.bef.TT0
beta.TT0 <- -B.TT0
beta.aft.TT0 <- -B.aft.TT0

# Calculate the mean elongation of the moon from the sun, X0.bef.TT0, X0.TT0, X0.aft.TT0 (in degrees)
X0.bef.TT0 <- 297.85036 + 445267.111480*JCE.bef.TT0 - 0.0019142*(JCE.bef.TT0^2) + (JCE.bef.TT0^3)/189474
X0.TT0 <- 297.85036 + 445267.111480*JCE.TT0 - 0.0019142*(JCE.TT0^2) + (JCE.TT0^3)/189474
X0.aft.TT0 <- 297.85036 + 445267.111480*JCE.aft.TT0 - 0.0019142*(JCE.aft.TT0^2) + (JCE.aft.TT0^3)/189474

# Calculate the mean anomaly of the sun (Earth), X1.bef.TT0, X1.TT0, X1.aft.TT0 (in degrees)
X1.bef.TT0 <- 357.52772 + 35999.050340*JCE.bef.TT0 - 0.0001603*(JCE.bef.TT0^2) - (JCE.bef.TT0^3)/300000
X1.TT0 <- 357.52772 + 35999.050340*JCE.TT0 - 0.0001603*(JCE.TT0^2) - (JCE.TT0^3)/300000
X1.aft.TT0 <- 357.52772 + 35999.050340*JCE.aft.TT0 - 0.0001603*(JCE.aft.TT0^2) - (JCE.aft.TT0^3)/300000

# Calculate the mean anomaly of the moon, X2.bef.TT0, X2.TT0, X2.aft.TT0 (in degrees)
X2.bef.TT0 <- 134.96298 + 477198.867398*JCE.bef.TT0 + 0.0086972*(JCE.bef.TT0^2) + (JCE.bef.TT0^3)/56250
X2.TT0 <- 134.96298 + 477198.867398*JCE.TT0 + 0.0086972*(JCE.TT0^2) + (JCE.TT0^3)/56250
X2.aft.TT0 <- 134.96298 + 477198.867398*JCE.aft.TT0 + 0.0086972*(JCE.aft.TT0^2) + (JCE.aft.TT0^3)/56250

# Calculate the moon's argument of latitude, X3.bef.TT0, X3.TT0, X3.aft.TT0 (in degrees)
X3.bef.TT0 <- 93.27191 + 483202.017538*JCE.bef.TT0 - 0.0036825*(JCE.bef.TT0^2) + (JCE.bef.TT0^3)/327270
X3.TT0 <- 93.27191 + 483202.017538*JCE.TT0 - 0.0036825*(JCE.TT0^2) + (JCE.TT0^3)/327270
X3.aft.TT0 <- 93.27191 + 483202.017538*JCE.aft.TT0 - 0.0036825*(JCE.aft.TT0^2) + (JCE.aft.TT0^3)/327270

# Calculate the longitude of the ascending node of the moon's mean orbit on the ecliptic,
# measured from the mean equinox of the date, X4.bef.TT0, X4.TT0, X4.aft.TT0 (in degrees)
X4.bef.TT0 <- 125.04452 - 1934.136261*JCE.bef.TT0 + 0.0020708*(JCE.bef.TT0^2) + (JCE.bef.TT0^3)/450000
X4.TT0 <- 125.04452 - 1934.136261*JCE.TT0 + 0.0020708*(JCE.TT0^2) + (JCE.TT0^3)/450000
X4.aft.TT0 <- 125.04452 - 1934.136261*JCE.aft.TT0 + 0.0020708*(JCE.aft.TT0^2) + (JCE.aft.TT0^3)/450000

# Calculate the terms Deltapsi.bef.TT0, Deltapsi.TT0, Deltapsi.aft.TT0 (in 0.0001 of arc seconds)
Deltapsi.bef.TT0 <- c()
Deltapsi.TT0 <- c()
Deltapsi.aft.TT0 <- c()
for (i in 1:dim(NLO)[1]){
  Deltapsi.bef.TT0[i] <- (NLO[i,6] + NLO[i,7]*JCE.bef.TT0) * sin(pi*(X0.bef.TT0*NLO[i,1]+X1.bef.TT0*NLO[i,2]+X2.bef.TT0*NLO[i,3]+X3.bef.TT0*NLO[i,4]+X4.bef.TT0*NLO[i,5])/180)
  Deltapsi.TT0[i] <- (NLO[i,6] + NLO[i,7]*JCE.TT0) * sin(pi*(X0.TT0*NLO[i,1]+X1.TT0*NLO[i,2]+X2.TT0*NLO[i,3]+X3.TT0*NLO[i,4]+X4.TT0*NLO[i,5])/180)
  Deltapsi.aft.TT0[i] <- (NLO[i,6] + NLO[i,7]*JCE.aft.TT0) * sin(pi*(X0.aft.TT0*NLO[i,1]+X1.aft.TT0*NLO[i,2]+X2.aft.TT0*NLO[i,3]+X3.aft.TT0*NLO[i,4]+X4.aft.TT0*NLO[i,5])/180)
}

# Calculate the nutation in longitude, Deltapsi.bef.TT0, Deltapsi.TT0, Deltapsi.aft.TT0 (in degrees)
Deltapsi.bef.TT0 <- sum(Deltapsi.bef.TT0)/36000000
Deltapsi.TT0 <- sum(Deltapsi.TT0)/36000000
Deltapsi.aft.TT0 <- sum(Deltapsi.aft.TT0)/36000000

# Calculate the terms Deltaepsilon.bef.TT0, Deltaepsilon.TT0, Deltaepsilon.aft.TT0 (in 0.0001 of arc seconds)
Deltaepsilon.bef.TT0 <- c()
Deltaepsilon.TT0 <- c()
Deltaepsilon.aft.TT0 <- c()
for (i in 1:dim(NLO)[1]){
  Deltaepsilon.bef.TT0[i] <- (NLO[i,8] + NLO[i,9]*JCE.bef.TT0) * cos(pi*(X0.bef.TT0*NLO[i,1]+X1.bef.TT0*NLO[i,2]+X2.bef.TT0*NLO[i,3]+X3.bef.TT0*NLO[i,4]+X4.bef.TT0*NLO[i,5])/180)
  Deltaepsilon.TT0[i] <- (NLO[i,8] + NLO[i,9]*JCE.TT0) * cos(pi*(X0.TT0*NLO[i,1]+X1.TT0*NLO[i,2]+X2.TT0*NLO[i,3]+X3.TT0*NLO[i,4]+X4.TT0*NLO[i,5])/180)
  Deltaepsilon.aft.TT0[i] <- (NLO[i,8] + NLO[i,9]*JCE.aft.TT0) * cos(pi*(X0.aft.TT0*NLO[i,1]+X1.aft.TT0*NLO[i,2]+X2.aft.TT0*NLO[i,3]+X3.aft.TT0*NLO[i,4]+X4.aft.TT0*NLO[i,5])/180)
}

# Calculate the nutation in obliquity, Deltaepsilon.bef.TT0, Deltaepsilon.TT0, Deltaepsilon.aft.TT0 (in degrees)
Deltaepsilon.bef.TT0 <- sum(Deltaepsilon.bef.TT0)/36000000
Deltaepsilon.TT0 <- sum(Deltaepsilon.TT0)/36000000
Deltaepsilon.aft.TT0 <- sum(Deltaepsilon.aft.TT0)/36000000

# Calculate the true obliquity of the ecliptic, epsilon0.bef.TT0, epsilon0.TT0, epsilon0.aft.TT0 (in arc seconds)
U.bef.TT0 <- JME.bef.TT0/10
epsilon0.bef.TT0 <- 84381.448 - 4680.93*U.bef.TT0 - 1.55*U.bef.TT0^2 + 1999.25*U.bef.TT0^3 - 51.38*U.bef.TT0^4 -249.67*U.bef.TT0^5 - 39.05*U.bef.TT0^6 + 7.12*U.bef.TT0^7 + 27.87*U.bef.TT0^8 + 5.79*U.bef.TT0^9 + 2.45*U.bef.TT0^10
U.TT0 <- JME.TT0/10
epsilon0.TT0 <- 84381.448 - 4680.93*U.TT0 - 1.55*U.TT0^2 + 1999.25*U.TT0^3 - 51.38*U.TT0^4 -249.67*U.TT0^5 - 39.05*U.TT0^6 + 7.12*U.TT0^7 + 27.87*U.TT0^8 + 5.79*U.TT0^9 + 2.45*U.TT0^10
U.aft.TT0 <- JME.aft.TT0/10
epsilon0.aft.TT0 <- 84381.448 - 4680.93*U.aft.TT0 - 1.55*U.aft.TT0^2 + 1999.25*U.aft.TT0^3 - 51.38*U.aft.TT0^4 -249.67*U.aft.TT0^5 - 39.05*U.aft.TT0^6 + 7.12*U.aft.TT0^7 + 27.87*U.aft.TT0^8 + 5.79*U.aft.TT0^9 + 2.45*U.aft.TT0^10

# Calculate the true obliquity of the ecliptic, epsilon.bef.TT0, epsilon.TT0, epsilon.aft.TT0 (in degrees) 
epsilon.bef.TT0 <- epsilon0.bef.TT0/3600 + Deltaepsilon.bef.TT0
epsilon.TT0 <- epsilon0.TT0/3600 + Deltaepsilon.TT0
epsilon.aft.TT0 <- epsilon0.aft.TT0/3600 + Deltaepsilon.aft.TT0

# Calculate the aberration correction, Deltatau.bef.TT0, Deltatau.TT0, Deltatau.aft.TT0 (in degrees):
Deltatau.bef.TT0 <- -20.4898/(3600*R.bef.TT0)
Deltatau.TT0 <- -20.4898/(3600*R.TT0)
Deltatau.aft.TT0 <- -20.4898/(3600*R.aft.TT0)

# Calculate the apparent sun longitude, lambda.bef.TT0, lambda.TT0, lambda.aft.TT0 (in degrees):
lambda.bef.TT0 <- Theta.limited.bef.TT0 + Deltapsi.bef.TT0 + Deltatau.bef.TT0
lambda.TT0 <- Theta.limited.TT0 + Deltapsi.TT0 + Deltatau.TT0
lambda.aft.TT0 <- Theta.limited.aft.TT0 + Deltapsi.aft.TT0 + Deltatau.aft.TT0


# Calculate the sun right ascension, alpha.bef.TT0, alpha.TT0, alpha.aft.TT0 (in radians)
alpha.bef.TT0 <- atan2((sin(lambda.bef.TT0*pi/180)*cos(epsilon.bef.TT0*pi/180)-tan(beta.bef.TT0*pi/180)*sin(epsilon.bef.TT0*pi/180)),cos(lambda.bef.TT0*pi/180))
alpha.TT0 <- atan2((sin(lambda.TT0*pi/180)*cos(epsilon.TT0*pi/180)-tan(beta.TT0*pi/180)*sin(epsilon.TT0*pi/180)),cos(lambda.TT0*pi/180))
alpha.aft.TT0 <- atan2((sin(lambda.aft.TT0*pi/180)*cos(epsilon.aft.TT0*pi/180)-tan(beta.aft.TT0*pi/180)*sin(epsilon.aft.TT0*pi/180)),cos(lambda.aft.TT0*pi/180))

# Calculate alpha.bef.TT0, alpha.TT0, alpha.aft.TT0 in degrees and limit it to the range from 0 to 360 degrees
alpha.bef.TT0 <- (alpha.bef.TT0*180)/pi
if (alpha.bef.TT0 < 0){
  alpha.limited.bef.TT0 <- 360 - 360*(abs(alpha.bef.TT0/360 - trunc(alpha.bef.TT0/360)))
}else{
  alpha.limited.bef.TT0 <- 360*(alpha.bef.TT0/360 - trunc(alpha.bef.TT0/360))
}
alpha.TT0 <- (alpha.TT0*180)/pi
if (alpha.TT0 < 0){
  alpha.limited.TT0 <- 360 - 360*(abs(alpha.TT0/360 - trunc(alpha.TT0/360)))
}else{
  alpha.limited.TT0 <- 360*(alpha.TT0/360 - trunc(alpha.TT0/360))
}
alpha.aft.TT0 <- (alpha.aft.TT0*180)/pi
if (alpha.aft.TT0 < 0){
  alpha.limited.aft.TT0 <- 360 - 360*(abs(alpha.aft.TT0/360 - trunc(alpha.aft.TT0/360)))
}else{
  alpha.limited.aft.TT0 <- 360*(alpha.aft.TT0/360 - trunc(alpha.aft.TT0/360))
}

# Calculate the geocentric sun declination, delta.bef.TT0, delta.TT0, delta.aft.TT0 (in degrees):
delta.bef.TT0 <- asin(sin(beta.bef.TT0*pi/180)*cos(epsilon.bef.TT0*pi/180)+cos(beta.bef.TT0*pi/180)*sin(epsilon.bef.TT0*pi/180)*sin(lambda.bef.TT0*pi/180))
delta.bef.TT0 <- (delta.bef.TT0*180)/pi
delta.TT0 <- asin(sin(beta.TT0*pi/180)*cos(epsilon.TT0*pi/180)+cos(beta.TT0*pi/180)*sin(epsilon.TT0*pi/180)*sin(lambda.TT0*pi/180))
delta.TT0 <- (delta.TT0*180)/pi
delta.aft.TT0 <- asin(sin(beta.aft.TT0*pi/180)*cos(epsilon.aft.TT0*pi/180)+cos(beta.aft.TT0*pi/180)*sin(epsilon.aft.TT0*pi/180)*sin(lambda.aft.TT0*pi/180))
delta.aft.TT0 <- (delta.aft.TT0*180)/pi


# Calculate the approximate sun transit time, m0, in fraction of day
m0 <- (alpha.limited.TT0 - sigma - v.UT0)/360

# Calculate the local hour angle corresponding to the sun elevation (h0Prime=-0.8333)
h0Prime <- -0.8333
arqument <- (sin(h0Prime*pi/180)-sin(phi*pi/180)*sin(delta.TT0*pi/180))/(cos(phi*pi/180)*cos(delta.TT0*pi/180))
if (abs(arqument) <=1){
  H0 <- acos(arqument)
}
else{
  warning('The sun is always above or below the horizon for that day')
}

H0 <- H0*180/pi

# Limit H0 to the range from 0 to 180 degrees
if (H0 < 0){
  H0.limited <- 180 - 180*(abs(H0/180 - trunc(H0/180)))
}else{
  H0.limited <- 180*(H0/180 - trunc(H0/180))
}

# Calculate the approximate sunset time, m1, in fraction of day
m1 <- m0 - H0.limited/360

# Calculate the approximate sunset time, m2, in fraction of day
m2 <- m0 + H0.limited/360

# Limit the value m0 to the range from 0 to 1 fraction of day
if (m0 < 0){
  m0.limited <- 1 - abs(m0 - trunc(m0))
}else{
  m0.limited <- m0 - trunc(m0)
}

# Limit the value m1 to the range from 0 to 1 fraction of day
if (m1 < 0){
  m1.limited <- 1 - abs(m1 - trunc(m1))
}else{
  m1.limited <- m1 - trunc(m1)
}

# Limit the value m2 to the range from 0 to 1 fraction of day
if (m2 < 0){
  m2.limited <- 1 - abs(m2 - trunc(m2))
}else{
  m2.limited <- m2 - trunc(m2)
}

# Calculate the sidereal time at Greenwich, in degrees, for the sun transit, sunrise and sunset, v0.TT0, v1.TT0, v2.TT0
v0.TT0 <- v.UT0 + 360.985647*m0.limited
v1.TT0 <- v.UT0 + 360.985647*m1.limited
v2.TT0 <- v.UT0 + 360.985647*m2.limited

# Calculate the terms n_i,
n0 <- m0.limited + DeltaT.TT0/86400
n1 <- m1.limited + DeltaT.TT0/86400
n2 <- m2.limited + DeltaT.TT0/86400

# Calculate the values alphaPrime_i and deltaPrime_i, in degrees
a <- alpha.limited.TT0-alpha.limited.bef.TT0
b <- alpha.limited.aft.TT0-alpha.limited.TT0
if (abs(a) > 2){
  if (a < 0){
    a <- 1 - abs(a - trunc(a))
  }else{
    a <- a - trunc(a)
  }
}
if (abs(b) > 2){
  if (b < 0){
    b <- 1 - abs(b - trunc(b))
  }else{
    b <- b - trunc(b)
  }
}
c <- b-a
aPrime <- delta.TT0-delta.bef.TT0
bPrime <- delta.aft.TT0-delta.TT0
if (abs(aPrime) > 2){
  if (aPrime < 0){
    aPrime <- 1 - abs(aPrime - trunc(aPrime))
  }else{
    aPrime <- aPrime - trunc(aPrime)
  }
}
if (abs(bPrime) > 2){
  if (bPrime < 0){
    bPrime <- 1 - abs(bPrime - trunc(bPrime))
  }else{
    bPrime <- bPrime - trunc(bPrime)
  }
}
cPrime <- bPrime-aPrime
alphaPrime0 <- alpha.limited.TT0 + (n0*(a+b+c*n0))/2
alphaPrime1 <- alpha.limited.TT0 + (n1*(a+b+c*n1))/2
alphaPrime2 <- alpha.limited.TT0 + (n2*(a+b+c*n2))/2
deltaPrime0 <- delta.TT0 + (n0*(aPrime+bPrime+cPrime*n0))/2
deltaPrime1 <- delta.TT0 + (n1*(aPrime+bPrime+cPrime*n1))/2
deltaPrime2 <- delta.TT0 + (n2*(aPrime+bPrime+cPrime*n2))/2

# Calculate the local hour angle for the sun transit, sunrise, and sunset, HPrime_i (in degrees)
HPrime0 <- v0.TT0 + sigma - alphaPrime0  
HPrime1 <- v1.TT0 + sigma - alphaPrime1  
HPrime2 <- v2.TT0 + sigma - alphaPrime2

# Limit HPrime0 to the range from -180 to 180 degrees
if (HPrime0 > 360){
  HPrime0.limited <- 360*(abs(HPrime0/360 - trunc(HPrime0/360)))
}else{
  if (HPrime0 < -360){
    HPrime0.limited <- -360*(abs(HPrime0/360 - trunc(HPrime0/360)))
  }else{
    HPrime0.limited <- HPrime0
  }	
}
HPrime0.limited
if (HPrime0.limited <= -180){
  HPrime0.limited <- HPrime0.limited + 360
}
if (HPrime0.limited >= 180){
  HPrime0.limited <- HPrime0.limited - 360
}

# Limit HPrime1 to the range from -180 to 180 degrees
if (HPrime1 > 360){
  HPrime1.limited <- 360*(abs(HPrime1/360 - trunc(HPrime1/360)))
}else{
  if (HPrime1 < -360){
    HPrime1.limited <- -360*(abs(HPrime1/360 - trunc(HPrime1/360)))
  }else{
    HPrime1.limited <- HPrime1
  }	
}
HPrime1.limited
if (HPrime1.limited <= -180){
  HPrime1.limited <- HPrime1.limited + 360
}
if (HPrime1.limited >= 180){
  HPrime1.limited <- HPrime1.limited - 360
}

# Limit HPrime2 to the range from -180 to 180 degrees
if (HPrime2 > 360){
  HPrime2.limited <- 360*(abs(HPrime2/360 - trunc(HPrime2/360)))
}else{
  if (HPrime2 < -360){
    HPrime2.limited <- -360*(abs(HPrime2/360 - trunc(HPrime2/360)))
  }else{
    HPrime2.limited <- HPrime2
  }	
}
HPrime2.limited
if (HPrime2.limited <= -180){
  HPrime2.limited <- HPrime2.limited + 360
}
if (HPrime2.limited >= 180){
  HPrime2.limited <- HPrime2.limited - 360
}

# Calculate the sun altitude for the sun transit, sunrise, and sunset, h_i (in degrees)
h0 <- asin(sin(phi*pi/180)*sin(deltaPrime0*pi/180)+cos(phi*pi/180)*cos(deltaPrime0*pi/180)*cos(HPrime0.limited*pi/180))
h1 <- asin(sin(phi*pi/180)*sin(deltaPrime1*pi/180)+cos(phi*pi/180)*cos(deltaPrime1*pi/180)*cos(HPrime1.limited*pi/180))
h2 <- asin(sin(phi*pi/180)*sin(deltaPrime2*pi/180)+cos(phi*pi/180)*cos(deltaPrime2*pi/180)*cos(HPrime2.limited*pi/180))

h0 <- h0*180/pi
h1 <- h1*180/pi
h2 <- h2*180/pi

# Calculate the sun transit, Transit (in fraction of day)
Transit <- m0.limited - HPrime0.limited/360

# Change to local time
Transit <- Transit + tz/24

# Limit Transit to the range from 0 to 1 
if (Transit < 0){
  Transit <- 1 - abs(Transit - trunc(Transit))
}else{
  Transit <- Transit - trunc(Transit)
}

# Change fraction of day to UT
Transit.hour <- trunc(24 * Transit)
Transit.minutes <- trunc(((24 * Transit) - trunc(24 * Transit)) * 60)
Transit.sec <- trunc((((24 * Transit) - trunc(24 * Transit)) * 60 - trunc(((24 * Transit) - trunc(24 * Transit)) * 60)) * 60)

Transit.time <- sprintf("%02d:%02d:%02d", Transit.hour, Transit.minutes, Transit.sec)

# Calculate the sunrise, Sunrise (in franction of day)
Sunrise <- m1.limited + (h1-h0Prime)/(360*cos(deltaPrime1*pi/180)*cos(phi*pi/180)*sin(HPrime1.limited*pi/180))

# Change to local time
Sunrise <- Sunrise + tz/24

# Limit Sunrise to the range from 0 to 1 
if (Sunrise < 0){
  Sunrise <- 1 - abs(Sunrise - trunc(Sunrise))
}else{
  Sunrise <- Sunrise - trunc(Sunrise)
}

# Change fraction of day to UT
Sunrise.hour <- trunc(24 * Sunrise)
Sunrise.minutes <- trunc(((24 * Sunrise) - trunc(24 * Sunrise)) * 60)
Sunrise.sec <- trunc((((24 * Sunrise) - trunc(24 * Sunrise)) * 60 - trunc(((24 * Sunrise) - trunc(24 * Sunrise)) * 60)) * 60)

Sunrise.time <- sprintf("%02d:%02d:%02d", Sunrise.hour, Sunrise.minutes, Sunrise.sec)

# Calculate the sunset, Sunset (in fraction of day)
Sunset <- m2.limited + (h2-h0Prime)/(360*cos(deltaPrime2*pi/180)*cos(phi*pi/180)*sin(HPrime2.limited*pi/180))

# Change to local time
Sunset <- Sunset + tz/24

# Limit Sunset to the range from 0 to 1 
if (Sunset < 0){
  Sunset <- 1 - abs(Sunset - trunc(Sunset))
}else{
  Sunset <- Sunset - trunc(Sunset)
}

# Change fraction of day to UT
Sunset.hour <- trunc(24 * Sunset)
Sunset.minutes <- trunc(((24 * Sunset) - trunc(24 * Sunset)) * 60)
Sunset.sec <- trunc((((24 * Sunset) - trunc(24 * Sunset)) * 60 - trunc(((24 * Sunset) - trunc(24 * Sunset)) * 60)) * 60)

Sunset.time <- sprintf("%02d:%02d:%02d", Sunset.hour, Sunset.minutes, Sunset.sec)

return(list(JD=JD,L=L.limited,B=B,R=R,Deltapsi=Deltapsi,Deltaepsilon=Deltaepsilon,epsilon=epsilon,H=H.limited,Zenith=theta,Azimuth=Phi.limited,Incidence=I,Sunrise.time=Sunrise.time, Transit.time=Transit.time, Sunset.time=Sunset.time))

}

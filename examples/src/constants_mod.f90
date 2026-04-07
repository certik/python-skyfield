! ─────────────────────────────────────────────────────────────────────
!  constants_mod — physical and mathematical constants
!  Superset of all program-specific constants_mod versions.
! ─────────────────────────────────────────────────────────────────────
module constants_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: PI       = 3.14159265358979323846264338327950288_dp
  real(dp), parameter :: TAU      = 2.0_dp * PI
  real(dp), parameter :: DEG2RAD  = PI / 180.0_dp
  real(dp), parameter :: RAD2DEG  = 180.0_dp / PI
  real(dp), parameter :: T0       = 2451545.0_dp           ! J2000.0 epoch (JD)
  real(dp), parameter :: DAY_S    = 86400.0_dp             ! seconds per day
  real(dp), parameter :: C_AUDAY  = 173.14463267424034_dp  ! speed of light AU/day
  real(dp), parameter :: C_SI     = 299792458.0_dp         ! speed of light m/s
  real(dp), parameter :: C_LIGHT_KMS = C_SI / 1000.0_dp   ! speed of light km/s
  real(dp), parameter :: AU_M     = 149597870700.0_dp      ! AU in metres
  real(dp), parameter :: AU_KM    = AU_M / 1000.0_dp       ! AU in km
  real(dp), parameter :: ERAD     = 6378136.6_dp           ! Earth equatorial radius (m)
  real(dp), parameter :: R_EARTH_KM = 6378.137_dp          ! Earth radius (km, WGS84 semi-major)
  real(dp), parameter :: ANGVEL   = 7.2921150d-5           ! Earth rotation rate rad/s
  real(dp), parameter :: ASEC2RAD = 4.848136811095359935899141d-6
  real(dp), parameter :: ASEC360  = 1296000.0_dp
  real(dp), parameter :: GS       = 1.32712440017987d+20   ! GM_sun m^3/s^2
  real(dp), parameter :: UNIX_JD_EPOCH = 2440587.5_dp      ! JD of Unix epoch (1970-01-01)

  ! J2000 obliquity for equatorial → ecliptic rotation
  real(dp), parameter :: OBLIQUITY_DEG = 23.4392794_dp
  real(dp), parameter :: OBLIQUITY_RAD = OBLIQUITY_DEG * PI / 180.0_dp

  ! WGS84
  real(dp), parameter :: WGS84_RADIUS = 6378137.0_dp
  real(dp), parameter :: WGS84_INVF   = 298.257223563_dp

  ! Reciprocal masses for deflection
  real(dp), parameter :: RMASS_SUN     = 1.0_dp
  real(dp), parameter :: RMASS_JUPITER = 1047.3486_dp
  real(dp), parameter :: RMASS_SATURN  = 3497.898_dp
  real(dp), parameter :: RMASS_EARTH   = 332946.050895_dp

  real(dp), parameter :: TENTH_USEC_2_RAD = ASEC2RAD / 1.0d7

  ! Body radii
  real(dp), parameter :: SOLAR_RADIUS_KM = 695700.0_dp
  real(dp), parameter :: MOON_RADIUS_KM  = 1737.4_dp
end module constants_mod

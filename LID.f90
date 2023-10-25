MODULE LID

   integer, parameter :: SURF = 1
   integer, parameter :: SOIL = 2
   integer, parameter :: STOR = 3


   ! LID Surface Layer
   type :: TSurfaceLayer
      real(kind=8) :: thickness     ! depression storage or berm ht. (ft)
      real(kind=8) :: voidFrac      ! available fraction of storage volume
      real(kind=8) :: roughness     ! surface Mannings n
      real(kind=8) :: surfSlope     ! land surface slope (fraction)
      real(kind=8) :: sideSlope     ! swale side slope (run/rise)
      real(kind=8) :: alpha         ! slope/roughness term in Manning eqn.
      logical :: canOverflow        ! .true. if immediate outflow of excess water
   end type TSurfaceLayer

! LID Pavement Layer
   type :: TPavementLayer
      real(kind=8) :: thickness     ! layer thickness (ft)
      real(kind=8) :: voidFrac      ! void volume / total volume
      real(kind=8) :: impervFrac    ! impervious area fraction
      real(kind=8) :: kSat          ! permeability (ft/sec)
      real(kind=8) :: clogFactor    ! clogging factor
      real(kind=8) :: regenDays     ! clogging regeneration interval (days)
      real(kind=8) :: regenDegree   ! degree of clogging regeneration
   end type TPavementLayer

! LID Soil Layer
   type :: TSoilLayer
      real(kind=8) :: thickness     ! layer thickness (ft)
      real(kind=8) :: porosity      ! void volume / total volume
      real(kind=8) :: fieldCap      ! field capacity
      real(kind=8) :: wiltPoint     ! wilting point
      real(kind=8) :: suction       ! suction head at wetting front (ft)
      real(kind=8) :: kSat          ! saturated hydraulic conductivity (ft/sec)
      real(kind=8) :: kSlope        ! slope of log(K) v. moisture content curve
   end type TSoilLayer

! LID Storage Layer
   type :: TStorageLayer
      real(kind=8) :: thickness     ! layer thickness (ft)
      real(kind=8) :: voidFrac      ! void volume / total volume
      real(kind=8) :: kSat          ! saturated hydraulic conductivity (ft/sec)
      real(kind=8) :: clogFactor    ! clogging factor
      logical :: covered            ! .true. if rain barrel is covered
   end type TStorageLayer

! Underdrain System (part of Storage Layer)
   type :: TDrainLayer
      real(kind=8) :: coeff         ! underdrain flow coeff. (in/hr or mm/hr)
      real(kind=8) :: expon         ! underdrain head exponent (for in or mm)
      real(kind=8) :: offset        ! offset height of underdrain (ft)
      real(kind=8) :: delay         ! rain barrel drain delay time (sec)
      real(kind=8) :: hOpen         ! head when drain opens (ft)
      real(kind=8) :: hClose        ! head when drain closes (ft)
      integer :: qCurve             ! curve controlling flow rate (optional)
   end type TDrainLayer

! Drainage Mat Layer (for green roofs)
   type :: TDrainMatLayer
      real(kind=8) :: thickness     ! layer thickness (ft)
      real(kind=8) :: voidFrac      ! void volume / total volume
      real(kind=8) :: roughness     ! Mannings n for green roof drainage mats
      real(kind=8) :: alpha         ! slope/roughness term in Manning equation
   end type TDrainMatLayer

! LID Process - generic LID design per unit of area
   type :: TLidProc
      character(len=:), allocatable :: ID    ! identifying name
      integer :: lidType                     ! type of LID
      type(TSurfaceLayer) :: surface         ! surface layer parameters
      type(TPavementLayer) :: pavement       ! pavement layer parameters
      type(TSoilLayer) :: soil               ! soil layer parameters
      type(TStorageLayer) :: storage         ! storage layer parameters
      type(TDrainLayer) :: drain             ! underdrain system parameters
      type(TDrainMatLayer) :: drainMat       ! drainage mat layer
      real(kind=8), allocatable :: drainRmvl(:) ! underdrain pollutant removals
   end type TLidProc

! Water Balance Statistics
   type :: TWaterBalance
      real(kind=8) :: inflow       ! total inflow (ft)
      real(kind=8) :: evap         ! total evaporation (ft)
      real(kind=8) :: infil        ! total infiltration (ft)
      real(kind=8) :: surfFlow     ! total surface runoff (ft)
      real(kind=8) :: drainFlow    ! total underdrain flow (ft)
      real(kind=8) :: initVol      ! initial stored volume (ft)
      real(kind=8) :: finalVol     ! final stored volume (ft)
   end type TWaterBalance

! LID Report File
   type :: TLidRptFile
      integer :: wasDry            ! number of successive dry periods
      character(len=256) :: results ! results for current time period
   end type TLidRptFile


! Green-Ampt Infiltration
   type :: TGrnAmpt
      real(kind=8) :: S               ! avg. capillary suction (ft)
      real(kind=8) :: Ks              ! saturated conductivity (ft/sec)
      real(kind=8) :: IMDmax          ! max. soil moisture deficit (ft/ft)
      !-----------------------------
      real(kind=8) :: IMD             ! current initial soil moisture deficit
      real(kind=8) :: F               ! current cumulative infiltrated volume (ft)
      real(kind=8) :: Fu              ! current upper zone infiltrated volume (ft)
      real(kind=8) :: Lu              ! depth of upper soil zone (ft)
      real(kind=8) :: T               ! time until start of next rain event (sec)
      logical :: Sat                  ! saturation flag
   end type TGrnAmpt

   
! LID Unit - specific LID process applied over a given area
   type :: TLidUnit
      integer :: lidIndex          ! index of LID process
      integer :: number            ! number of replicate units
      real(kind=8) :: area         ! area of single replicate unit (ft2)
      real(kind=8) :: fullWidth    ! full top width of single unit (ft)
      real(kind=8) :: botWidth     ! bottom width of single unit (ft)
      real(kind=8) :: initSat      ! initial saturation of soil & storage layers
      real(kind=8) :: fromImperv   ! fraction of impervious area runoff treated
      real(kind=8) :: fromPerv     ! fraction of pervious area runoff treated
      logical :: toPerv            ! .true. if outflow sent to pervious area; .false. if not
      integer :: drainSubcatch     ! subcatchment receiving drain flow
      integer :: drainNode         ! node receiving drain flow
      type(TLidRptFile) :: rptFile ! pointer to detailed report file

      type(TGrnAmpt) :: soilInfil  ! infil. object for biocell soil layer
      real(kind=8) :: surfaceDepth ! depth of ponded water on surface layer (ft)
      real(kind=8) :: paveDepth    ! depth of water in porous pavement layer
      real(kind=8) :: soilMoisture ! moisture content of biocell soil layer
      real(kind=8) :: storageDepth ! depth of water in storage layer (ft)

! net inflow - outflow from previous time step for each LID layer (ft/s)
      real(kind=8) :: oldFluxRates(MAX_LAYERS)

      real(kind=8) :: dryTime       ! time since last rainfall (sec)
      real(kind=8) :: oldDrainFlow  ! previous drain flow (cfs)
      real(kind=8) :: newDrainFlow
   end type TLidUnit



contains


   subroutine greenRoofFluxRates(x, f)
      implicit none

      ! 输入: x = 存储水平的向量
      ! 输出: f = 通量速率的向量

      real(kind=8), intent(in) :: x(:)
      real(kind=8), intent(out) :: f(:)

      ! 湿度水平变量
      real(kind=8) :: surfaceDepth
      real(kind=8) :: soilTheta
      real(kind=8) :: storageDepth

      ! 中间变量
      real(kind=8) :: availVolume
      real(kind=8) :: maxRate

      ! 绿色屋顶属性
      real(kind=8) :: soilThickness
      real(kind=8) :: storageThickness
      real(kind=8) :: soilPorosity
      real(kind=8) :: storageVoidFrac
      real(kind=8) :: soilFieldCap
      real(kind=8) :: soilWiltPoint

      !... 从输入向量中获取湿度水平
      surfaceDepth = x(SURF)
      soilTheta = x(SOIL)
      storageDepth = x(STOR)

      !... 将湿度水平转换为体积
      SurfaceVolume = surfaceDepth * theLidProc%surface%voidFrac
      SoilVolume = soilTheta * soilThickness
      StorageVolume = storageDepth * storageVoidFrac

      !... 获取蒸散发速率
      availVolume = SoilVolume - soilWiltPoint * soilThickness
      call getEvapRates(SurfaceVolume, 0.0, availVolume, StorageVolume, 1.0)
      if (soilTheta >= soilPorosity) StorageEvap = 0.0

      !... 土壤层渗透速率
      SoilPerc = getSoilPercRate(soilTheta)

      !... 将渗透速率限制为可用水量
      availVolume = (soilTheta - soilFieldCap) * soilThickness
      maxRate = MAX(availVolume, 0.0) / Tstep - SoilEvap
      SoilPerc = MIN(SoilPerc, maxRate)
      SoilPerc = MAX(SoilPerc, 0.0)

      !... 蓄水层（排水层）流出速率
      StorageExfil = 0.0
      StorageDrain = getDrainMatOutflow(storageDepth)

      !... 单元已满
      if (soilTheta >= soilPorosity .and. storageDepth >= storageThickness) then
         !... 两层的流出速率等于限制速率
         maxRate = MIN(SoilPerc, StorageDrain)
         SoilPerc = maxRate
         StorageDrain = maxRate

         !... 调整入流速率到土壤层
         SurfaceInfil = MIN(SurfaceInfil, maxRate)
      else
         !... 限制排水层流出速率为可用储存体积
         maxRate = storageDepth * storageVoidFrac / Tstep - StorageEvap
         if (storageDepth >= storageThickness) maxRate = maxRate + SoilPerc
         maxRate = MAX(maxRate, 0.0)
         StorageDrain = MIN(StorageDrain, maxRate)

         !... 限制土壤渗透入流速率为未使用的储存体积
         maxRate = (storageThickness - storageDepth) * storageVoidFrac / Tstep + &
            StorageDrain + StorageEvap
         SoilPerc = MIN(SoilPerc, maxRate)

         !... 调整地表入渗速率，使得土壤孔隙度不超过上限
         maxRate = (soilPorosity - soilTheta) * soilThickness / Tstep + &
            SoilPerc + SoilEvap
         SurfaceInfil = MIN(SurfaceInfil, maxRate)
      end if

      ! ... 计算地表流出速率
      SurfaceOutflow = getSurfaceOutflowRate(surfaceDepth)

      ! ... 计算整体层的通量速率
      f(SURF) = (SurfaceInflow - SurfaceEvap - SurfaceInfil - SurfaceOutflow) / &
         theLidProc%surface%voidFrac
      f(SOIL) = (SurfaceInfil - SoilEvap - SoilPerc) / &
         theLidProc%soil%thickness
      f(STOR) = (SoilPerc - StorageEvap - StorageDrain) / &
         theLidProc%storage%voidFrac


   end

END MODULE LID

! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_paraview_catalyst

   use mpas_grid_types
   use mpas_dmpar
   use mpas_dmpar_types
   use mpas_sort
   use mpas_configure
   use mpas_io_streams
   use mpas_io_output

   contains

   subroutine mpas_insitu_initialize!{{{

      call mpas_initialize

   end subroutine mpas_insitu_initialize!}}}

   subroutine mpas_insitu_create_geometry(domain)!{{{
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (block_type), pointer :: block_ptr
      type (mpas_exchange_list), pointer :: exchListPtr
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      integer :: nCells, nEdges, nVertices, maxEdges, vertexDegree, nVertLevels
      integer :: nGhostCell, nGhostVertex
      integer :: nHaloLayers, iHalo, iCell, indx
      integer :: indexT, indexS, numTracers
      integer, dimension(:), allocatable :: cellGhost, cellHalo
      integer, dimension(:), allocatable :: vertexGhost, vertexHalo

      block_ptr => domain % blocklist

      nCells = block_ptr % mesh % nCells
      nVertices = block_ptr % mesh % nVertices
      nVertLevels = block_ptr % mesh % nVertLevels
      nEdges = block_ptr % mesh % nEdges
      maxEdges = block_ptr % mesh % maxEdges
      vertexDegree = block_ptr % mesh % vertexDegree

      ! Collect the number of halo cells in all levels
      nHaloLayers = config_num_halos
      nGhostCell = 0
      do iHalo = 1, nHaloLayers
        exchListPtr => block_ptr %parinfo % cellsToRecv % halos(iHalo) % exchList
        do while(associated(exchListPtr))
          nGhostCell = nGhostCell + exchListPtr % nList
          exchListPtr => exchListPtr % next
        end do
      end do

      nGhostVertex = 0
      do iHalo = 1, nHaloLayers
        exchListPtr => block_ptr %parinfo % verticesToRecv % halos(iHalo) % exchList
        do while(associated(exchListPtr))
          nGhostVertex = nGhostVertex + exchListPtr % nList
          exchListPtr => exchListPtr % next
        end do
      end do

      ! Allocate arrays to hold halo cell index and halo number
      allocate(cellGhost(nGhostCell))
      allocate(cellHalo(nGhostCell))
      allocate(vertexGhost(nGhostVertex))
      allocate(vertexHalo(nGhostVertex))

      ! Fill in halo cells indices and halo number
      indx = 1
      do iHalo = 1, nHaloLayers
        exchListPtr => block_ptr %parinfo % cellsToRecv % halos(iHalo) % exchList
        do while(associated(exchListPtr))
          do iCell = 1, exchListPtr % nList
            cellGhost(indx) = exchListPtr % destList(iCell)
            cellHalo(indx) = iHalo
            indx = indx + 1
          end do
          exchListPtr => exchListPtr % next
        end do
      end do

      indx = 1
      do iHalo = 1, nHaloLayers
        exchListPtr => block_ptr %parinfo % verticesToRecv % halos(iHalo) % exchList
        do while(associated(exchListPtr))
          do iCell = 1, exchListPtr % nList
            vertexGhost(indx) = exchListPtr % destList(iCell)
            vertexHalo(indx) = iHalo
            indx = indx + 1
          end do
          exchListPtr => exchListPtr % next
        end do
      end do

      print *,'COPROCESS create geometry',nCells,nVertices,maxEdges,vertexDegree,nVertLevels
      call coprocessor_create_grid( \
                       nCells, maxEdges, nGhostCell, cellGhost, cellHalo, \
                       nVertices, vertexDegree, nGhostVertex, vertexGhost, vertexHalo, \
                       nVertLevels, \
                       block_ptr % mesh % xCell % array, \
                       block_ptr % mesh % yCell % array, \
                       block_ptr % mesh % zCell % array, \
                       block_ptr % mesh % xVertex % array, \
                       block_ptr % mesh % yVertex % array, \
                       block_ptr % mesh % zVertex % array, \
                       block_ptr % mesh % lonCell % array, \
                       block_ptr % mesh % latCell % array, \
                       block_ptr % mesh % lonVertex % array, \
                       block_ptr % mesh % latVertex % array, \
                       block_ptr % mesh % nEdgesOnCell % array, \
                       block_ptr % mesh % cellsOnVertex % array, \
                       block_ptr % mesh % vertexMask % array, \
                       block_ptr % mesh % verticesOnCell % array, \
                       block_ptr % mesh % cellMask % array)

      tracers => block_ptr % state % time_levs(2) % state % tracers % array
      numTracers = size(tracers, dim=1)
      indexT = block_ptr % state % time_levs(1) % state % index_temperature
      indexS = block_ptr % state % time_levs(1) % state % index_salinity
      call coprocessor_register_tracer_data( \
                       indexT, \
                       block_ptr % state % time_levs(1) % state % tracers % constituentNames(indexT), \
                       numTracers, \
                       block_ptr % state % time_levs(1) % state % tracers % dimSizes(2), \
                       block_ptr % state % time_levs(1) % state % tracers % dimSizes(3), \
                       block_ptr % state % time_levs(1) % state % tracers % array)
      call coprocessor_register_tracer_data( \
                       indexS, \
                       block_ptr % state % time_levs(1) % state % tracers % constituentNames(indexS), \
                       numTracers, \
                       block_ptr % state % time_levs(1) % state % tracers % dimSizes(2), \
                       block_ptr % state % time_levs(1) % state % tracers % dimSizes(3), \
                       block_ptr % state % time_levs(1) % state % tracers % array)

      call coprocessor_register_data( \
                       block_ptr % diagnostics % density % fieldName, \
                       block_ptr % diagnostics % density % dimSizes(1), \
                       block_ptr % diagnostics % density % dimSizes(2), \
                       block_ptr % diagnostics % density % array)
      call coprocessor_register_data( \
                       block_ptr % diagnostics % pressure % fieldName, \
                       block_ptr % diagnostics % pressure % dimSizes(1), \
                       block_ptr % diagnostics % pressure % dimSizes(2), \
                       block_ptr % diagnostics % pressure % array)
      call coprocessor_register_data( \
                       block_ptr % diagnostics % circulation % fieldName, \
                       block_ptr % diagnostics % circulation % dimSizes(1), \
                       block_ptr % diagnostics % circulation % dimSizes(2), \
                       block_ptr % diagnostics % circulation % array)
      call coprocessor_register_data( \
                       block_ptr % diagnostics % relativeVorticity % fieldName, \
                       block_ptr % diagnostics % relativeVorticity % dimSizes(1), \
                       block_ptr % diagnostics % relativeVorticity % dimSizes(2), \
                       block_ptr % diagnostics % relativeVorticity % array)

      deallocate(cellGhost)
      deallocate(cellHalo)
      deallocate(vertexGhost)
      deallocate(vertexHalo)

   end subroutine mpas_insitu_create_geometry!}}}

   subroutine mpas_insitu_load_data(domain, itime)!{{{

      implicit none
 
      type (domain_type), intent(inout) :: domain
      type (block_type), pointer :: block_ptr
      integer, intent(inout) :: itime
      integer :: indexT, indexS, numTracers
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      block_ptr => domain % blocklist

      tracers => block_ptr % state % time_levs(2) % state % tracers % array
      numTracers = size(tracers, dim=1)
      indexT = block_ptr % state % time_levs(1) % state % index_temperature
      indexS = block_ptr % state % time_levs(1) % state % index_salinity

      call coprocessor_add_tracer_data( \
                       itime, \
                       indexT, \
                       block_ptr % state % time_levs(1) % state % tracers % constituentNames(indexT), \
                       numTracers, \
                       block_ptr % state % time_levs(1) % state % tracers % dimSizes(2), \
                       block_ptr % state % time_levs(1) % state % tracers % dimSizes(3), \
                       block_ptr % state % time_levs(1) % state % tracers % array)
      call coprocessor_add_tracer_data( \
                       itime, \
                       indexS, \
                       block_ptr % state % time_levs(1) % state % tracers % constituentNames(indexS), \
                       numTracers, \
                       block_ptr % state % time_levs(1) % state % tracers % dimSizes(2), \
                       block_ptr % state % time_levs(1) % state % tracers % dimSizes(3), \
                       block_ptr % state % time_levs(1) % state % tracers % array)

      call coprocessor_add_data( \
                       itime, \
                       block_ptr % diagnostics % density % fieldName, \
                       block_ptr % diagnostics % density % dimSizes(1), \
                       block_ptr % diagnostics % density % dimSizes(2), \
                       block_ptr % diagnostics % density % array)
      call coprocessor_add_data( \
                       itime, \
                       block_ptr % diagnostics % pressure % fieldName, \
                       block_ptr % diagnostics % pressure % dimSizes(1), \
                       block_ptr % diagnostics % pressure % dimSizes(2), \
                       block_ptr % diagnostics % pressure % array)
      call coprocessor_add_data( \
                       itime, \
                       block_ptr % diagnostics % circulation % fieldName, \
                       block_ptr % diagnostics % circulation % dimSizes(1), \
                       block_ptr % diagnostics % circulation % dimSizes(2), \
                       block_ptr % diagnostics % circulation % array)
      call coprocessor_add_data( \
                       itime, \
                       block_ptr % diagnostics % relativeVorticity % fieldName, \
                       block_ptr % diagnostics % relativeVorticity % dimSizes(1), \
                       block_ptr % diagnostics % relativeVorticity % dimSizes(2), \
                       block_ptr % diagnostics % relativeVorticity % array)

   end subroutine mpas_insitu_load_data!}}}

   subroutine mpas_insitu_coprocess(domain, itime)!{{{

      implicit none

      type (domain_type), intent(inout) :: domain
      integer, intent(inout) :: itime
      integer :: doWork

      doWork = 0
      print *,'checking for COPROCESS ',itime
      call mpas_check_coprocess(itime, doWork)
      if(doWork .eq. 0) then
         return
      endif
      ! we actually need to do coprocessing
      ! the grids should already be loaded.
      call mpas_insitu_load_data(domain, itime)
      print *,'COPROCESS ',itime
      call mpas_coprocess

   end subroutine mpas_insitu_coprocess!}}}


   subroutine mpas_insitu_finalize()!{{{
 
      implicit none
 
      print *,'COPROCESS finalize'
      call mpas_finalize
 
   end subroutine mpas_insitu_finalize!}}}
 
end module mpas_paraview_catalyst

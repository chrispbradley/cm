!> \file
!> \author Chris Bradley
!> \brief This module contains all interface routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s): David Nordsletten, Thiranja Prasad Babarenda Gamage
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delte
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module contains all interface routines.
MODULE INTERFACE_ROUTINES

  USE BASE_ROUTINES
  USE FIELD_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE INTERFACE_LABEL_GET
    MODULE PROCEDURE INTERFACE_LABEL_GET_C
    MODULE PROCEDURE INTERFACE_LABEL_GET_VS
  END INTERFACE !INTERFACE_LABEL_GET
  
  INTERFACE INTERFACE_LABEL_SET
    MODULE PROCEDURE INTERFACE_LABEL_SET_C
    MODULE PROCEDURE INTERFACE_LABEL_SET_VS
  END INTERFACE !INTERFACE_LABEL_SET
  
  PUBLIC INTERFACE_MESH_ADD

  PUBLIC INTERFACE_CREATE_START, INTERFACE_CREATE_FINISH
  
  PUBLIC INTERFACE_SELF_CONTACT_SET
  
  PUBLIC INTERFACE_COORDINATE_SYSTEM_SET,INTERFACE_COORDINATE_SYSTEM_GET

  PUBLIC INTERFACE_DESTROY, INTERFACE_MESH_CONNECTIVITY_DESTROY
  
  PUBLIC INTERFACE_DATA_POINTS_GET

  PUBLIC INTERFACE_LABEL_GET,INTERFACE_LABEL_SET

  PUBLIC INTERFACE_MESH_CONNECTIVITY_CREATE_START, INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH

  PUBLIC INTERFACE_USER_NUMBER_FIND

  PUBLIC INTERFACES_FINALISE,INTERFACES_INITIALISE

  PUBLIC INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET, INTERFACE_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET

  PUBLIC INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_CONTACT_SET 
 
  PUBLIC INTERFACE_MESH_CONNECTIVITY_SET_BASIS
  
  PUBLIC INTERFACE_POINTS_CONNECTIVITY_CREATE_START, INTERFACE_POINTS_CONNECTIVITY_CREATE_FINISH
  
  PUBLIC INTERFACE_POINTS_CONNECTIVITY_ELEMENT_NUMBER_SET, INTERFACE_POINTS_CONNECTIVITY_POINT_XI_SET
  
  PUBLIC INTERFACE_POINTS_CONNECTIVITY_POINT_XI_CONTACT_SET
  
  PUBLIC INTERFACE_POINTS_CONNECTIVITY_PROJECTION_RESULTS_SET
  
CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE INTERFACE_MESH_ADD(INTERFACE,MESH,MESH_INDEX,ERR,ERROR,*)   

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to add a mesh to
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to add to the interface
    INTEGER(INTG), INTENT(OUT) :: MESH_INDEX !<On return, the index of the added mesh in the list of meshes in the interface
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx
    LOGICAL :: MESH_ALREADY_COUPLED
    TYPE(MESH_TYPE), POINTER :: COUPLED_MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_COUPLED_MESHES(:)
    TYPE(REGION_TYPE), POINTER :: COUPLED_MESH_REGION,MESH_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_COUPLED_MESHES)
    
    CALL ENTERS("INTERFACE_MESH_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FLAG_ERROR("Interface has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(MESH)) THEN
          IF(MESH%MESH_FINISHED) THEN
            MESH_REGION=>MESH%REGION
            IF(ASSOCIATED(MESH_REGION)) THEN
              ALLOCATE(NEW_COUPLED_MESHES(INTERFACE%NUMBER_OF_COUPLED_MESHES+1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new coupled meshes.",ERR,ERROR,*999)
              !Check that the mesh is not already in the list of meshes for the interface.
              IF(INTERFACE%NUMBER_OF_COUPLED_MESHES>0) THEN
                IF(ASSOCIATED(INTERFACE%COUPLED_MESHES)) THEN
                  MESH_ALREADY_COUPLED=.FALSE.
                  DO mesh_idx=1,INTERFACE%NUMBER_OF_COUPLED_MESHES
                    COUPLED_MESH=>INTERFACE%COUPLED_MESHES(mesh_idx)%PTR
                    IF(ASSOCIATED(COUPLED_MESH)) THEN
                      COUPLED_MESH_REGION=>COUPLED_MESH%REGION
                      IF(ASSOCIATED(COUPLED_MESH_REGION)) THEN
                        IF(MESH_REGION%USER_NUMBER==COUPLED_MESH_REGION%USER_NUMBER) THEN
                          IF(MESH%USER_NUMBER==COUPLED_MESH%USER_NUMBER .AND. .NOT. INTERFACE%SELF_CONTACT) THEN
                            MESH_ALREADY_COUPLED=.TRUE.
                            EXIT
                          ENDIF
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Coupled interface mesh region for mesh index "// &
                          & TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))//" is not associated."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Coupled interface mesh for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                        & " is not associated."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    NEW_COUPLED_MESHES(mesh_idx)%PTR=>INTERFACE%COUPLED_MESHES(mesh_idx)%PTR
                  ENDDO !mesh_idx
                  IF(MESH_ALREADY_COUPLED) THEN
                    LOCAL_ERROR="The supplied mesh has already been added to the list of coupled meshes at mesh index "// &
                      & TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  DEALLOCATE(INTERFACE%COUPLED_MESHES)
                ELSE
                  CALL FLAG_ERROR("Interface coupled meshes is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              !Add the mesh to the list of coupled meshes
              NEW_COUPLED_MESHES(INTERFACE%NUMBER_OF_COUPLED_MESHES+1)%PTR=>MESH
              INTERFACE%COUPLED_MESHES=>NEW_COUPLED_MESHES
              !Increment the number of coupled meshes and return the index
              INTERFACE%NUMBER_OF_COUPLED_MESHES=INTERFACE%NUMBER_OF_COUPLED_MESHES+1
              MESH_INDEX=INTERFACE%NUMBER_OF_COUPLED_MESHES
            ELSE
              CALL FLAG_ERROR("Mesh region is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MESH_ADD")
    RETURN
999 IF(ASSOCIATED(NEW_COUPLED_MESHES)) DEALLOCATE(NEW_COUPLED_MESHES)
    CALL ERRORS("INTERFACE_MESH_ADD",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_ADD")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_ADD

  !
  !================================================================================================================================
  !

  !>Finishes the creation of an interface. \see OPENCMISS::CMISSInterfaceCreateFinish
  SUBROUTINE INTERFACE_CREATE_FINISH(INTERFACE,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    CALL ENTERS("INTERFACE_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FLAG_ERROR("Interface has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(.NOT. INTERFACE%SELF_CONTACT) THEN
          IF(INTERFACE%NUMBER_OF_COUPLED_MESHES<2) THEN
            LOCAL_ERROR="Invalid mesh coupling. Only "//TRIM(NUMBER_TO_VSTRING(INTERFACE%NUMBER_OF_COUPLED_MESHES, &
              & "*",ERR,ERROR))// " have been coupled. The number of coupled meshes must be >= 2."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
        INTERFACE%INTERFACE_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Interface :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  User number = ",INTERFACE%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Global number = ",INTERFACE%GLOBAL_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",INTERFACE%LABEL,ERR,ERROR,*999)
      IF(ASSOCIATED(INTERFACE%INTERFACES)) THEN
        IF(ASSOCIATED(INTERFACE%INTERFACES%PARENT_REGION)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region user number = ",INTERFACE%INTERFACES% &
            & PARENT_REGION%USER_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region label = ",INTERFACE%INTERFACES% &
            & PARENT_REGION%LABEL,ERR,ERROR,*999)        
        ELSE
          CALL FLAG_ERROR("Interfaces parent region is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface interfaces is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    CALL EXITS("INTERFACE_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE INTERFACE_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of an interface on a parent region. \see OPENCMISS::CMISSInterfaceCreateStart
  SUBROUTINE INTERFACE_CREATE_START(USER_NUMBER,PARENT_REGION,INTERFACE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the interface to create
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<A pointer to the parent region to create the interface on.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<On exit, a pointer to the created interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_idx
    TYPE(INTERFACE_TYPE), POINTER :: NEW_INTERFACE
    TYPE(INTERFACE_PTR_TYPE), POINTER :: NEW_INTERFACES(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR,LOCAL_STRING

    NULLIFY(NEW_INTERFACE)
    NULLIFY(NEW_INTERFACES)
    
    CALL ENTERS("INTERFACE_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(PARENT_REGION)) THEN
      IF(ASSOCIATED(INTERFACE)) THEN
        CALL FLAG_ERROR("Interface is already associated.",ERR,ERROR,*998)
      ELSE
        NULLIFY(INTERFACE)
        CALL INTERFACE_USER_NUMBER_FIND(USER_NUMBER,PARENT_REGION,INTERFACE,ERR,ERROR,*998)
        IF(ASSOCIATED(INTERFACE)) THEN
          LOCAL_ERROR="Interface number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
        ELSE        
          NULLIFY(INTERFACE)
          !Allocate and set default interface properties.
          CALL INTERFACE_INITIALISE(NEW_INTERFACE,ERR,ERROR,*999)
          NEW_INTERFACE%USER_NUMBER=USER_NUMBER
          NEW_INTERFACE%GLOBAL_NUMBER=PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES+1
          LOCAL_STRING="Interface_"//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
          NEW_INTERFACE%LABEL=CHAR(LOCAL_STRING)
          IF(ERR/=0) GOTO 999
          NEW_INTERFACE%INTERFACES=>PARENT_REGION%INTERFACES
          NEW_INTERFACE%PARENT_REGION=>PARENT_REGION
          NEW_INTERFACE%SELF_CONTACT=.FALSE.
          !Add new initerface into list of interfaces in the parent region
          ALLOCATE(NEW_INTERFACES(PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES+1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new interfaces.",ERR,ERROR,*999)
          DO interface_idx=1,PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES
            NEW_INTERFACES(interface_idx)%PTR=>PARENT_REGION%INTERFACES%INTERFACES(interface_idx)%PTR
          ENDDO !interface_idx
          NEW_INTERFACES(PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES+1)%PTR=>NEW_INTERFACE
          IF(ASSOCIATED(PARENT_REGION%INTERFACES%INTERFACES)) DEALLOCATE(PARENT_REGION%INTERFACES%INTERFACES)
          PARENT_REGION%INTERFACES%INTERFACES=>NEW_INTERFACES
          PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES=PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES+1
          INTERFACE=>NEW_INTERFACE
          
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Parent region is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("INTERFACE_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_INTERFACES)) DEALLOCATE(NEW_INTERFACES)
    CALL INTERFACE_FINALISE(INTERFACE,ERR,ERROR,*998)
998 CALL ERRORS("INTERFACE_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_CREATE_START")
    RETURN 1
  END SUBROUTINE INTERFACE_CREATE_START

  !
  !================================================================================================================================
  !
  
  !>Finalises an interface and deallocates all memory.
  SUBROUTINE INTERFACE_SELF_CONTACT_SET(INTERFACE,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("INTERFACE_SELF_CONTACT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      INTERFACE%SELF_CONTACT=.TRUE.
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
      
    ENDIF
    
    CALL EXITS("INTERFACE_SELF_CONTACT_SET")
    RETURN
999 CALL ERRORS("INTERFACE_SELF_CONTACT_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_SELF_CONTACT_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_SELF_CONTACT_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the coordinate system of an interface. \see OPENCMISS::CMISSInterface_CoordinateSystemGet
  SUBROUTINE INTERFACE_COORDINATE_SYSTEM_GET(INTERFACE,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<On exit, the coordinate system for the specified interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_COORDINATE_SYSTEM_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          CALL FLAG_ERROR("Coordinate system is already associated.",ERR,ERROR,*999)
        ELSE
          COORDINATE_SYSTEM=>INTERFACE%COORDINATE_SYSTEM
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_COORDINATE_SYSTEM_GET")
    RETURN
999 CALL ERRORS("INTERFACE_COORDINATE_SYSTEM_GET",ERR,ERROR)
    CALL EXITS("INTERFACE_COORDINATE_SYSTEM_GET")
    RETURN 1
  END SUBROUTINE INTERFACE_COORDINATE_SYSTEM_GET
  
  
  !
  !================================================================================================================================
  !

  !>Sets the coordinate system of an interface.  \see OPENCMISS::CMISSInterface_CoordinateSystemSet
  SUBROUTINE INTERFACE_COORDINATE_SYSTEM_SET(INTERFACE,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to set the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<The coordinate system to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_COORDINATE_SYSTEM_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FLAG_ERROR("Interface has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
            INTERFACE%COORDINATE_SYSTEM=>COORDINATE_SYSTEM
          ELSE
            CALL FLAG_ERROR("Coordinate system has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_COORDINATE_SYSTEM_SET")
    RETURN
999 CALL ERRORS("INTERFACE_COORDINATE_SYSTEM_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_COORDINATE_SYSTEM_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_COORDINATE_SYSTEM_SET

  !
  !================================================================================================================================
  !

  !>Destroys an interface. \see OPENCMISS::CMISSInterfaceDestroy
  SUBROUTINE INTERFACE_DESTROY(INTERFACE,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_idx,interface_position
    TYPE(INTERFACE_PTR_TYPE), POINTER :: NEW_INTERFACES(:)
    TYPE(INTERFACES_TYPE), POINTER :: INTERFACES
     
    NULLIFY(NEW_INTERFACES)

    CALL ENTERS("INTERFACE_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      INTERFACES=>INTERFACE%INTERFACES
      IF(ASSOCIATED(INTERFACES)) THEN
        interface_position=INTERFACE%GLOBAL_NUMBER

        !Destroy all the interface condition components
        CALL INTERFACE_FINALISE(INTERFACE,ERR,ERROR,*999)
        
        !Remove the interface condition from the list of interface conditions
        IF(INTERFACES%NUMBER_OF_INTERFACES>1) THEN
          ALLOCATE(NEW_INTERFACES(INTERFACES%NUMBER_OF_INTERFACES-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new interface conditions.",ERR,ERROR,*999)
          DO interface_idx=1,INTERFACES%NUMBER_OF_INTERFACES
            IF(interface_idx<interface_position) THEN
              NEW_INTERFACES(interface_idx)%PTR=>INTERFACES%INTERFACES(interface_idx)%PTR
            ELSE IF(interface_idx>interface_position) THEN
              INTERFACES%INTERFACES(interface_idx)%PTR%GLOBAL_NUMBER=INTERFACES%INTERFACES(interface_idx)%PTR%GLOBAL_NUMBER-1
              NEW_INTERFACES(interface_idx-1)%PTR=>INTERFACES%INTERFACES(interface_idx)%PTR
            ENDIF
          ENDDO !interface_idx
          IF(ASSOCIATED(INTERFACES%INTERFACES)) DEALLOCATE(INTERFACES%INTERFACES)
          INTERFACES%INTERFACES=>NEW_INTERFACES
          INTERFACES%NUMBER_OF_INTERFACES=INTERFACES%NUMBER_OF_INTERFACES-1
        ELSE
          DEALLOCATE(INTERFACES%INTERFACES)
          INTERFACES%NUMBER_OF_INTERFACES=0
        ENDIF
        
      ELSE
        CALL FLAG_ERROR("Interface interfaces is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("INTERFACE_DESTROY")
    RETURN
999 CALL ERRORS("INTERFACE_DESTROY",ERR,ERROR)
    CALL EXITS("INTERFACE_DESTROY")
    RETURN 1
  END SUBROUTINE INTERFACE_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data points for a region. \see OPENCMISS::CMISSRegionDataPointsGet
  SUBROUTINE INTERFACE_DATA_POINTS_GET(INTERFACE,DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the region to get the data points for
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<On exit, a pointer to the data points for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("REGION_DATA_POINTS_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN 
        IF(ASSOCIATED(DATA_POINTS)) THEN
          CALL FLAG_ERROR("Data points is already associated.",ERR,ERROR,*998)
        ELSE
          DATA_POINTS=>INTERFACE%DATA_POINTS
          IF(.NOT.ASSOCIATED(DATA_POINTS)) CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_DATA_POINTS_GET")
    RETURN
999 NULLIFY(DATA_POINTS)
998 CALL ERRORS("INTERFACE_DATA_POINTS_GET",ERR,ERROR)
    CALL EXITS("INTERFACE_DATA_POINTS_GET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_DATA_POINTS_GET  

  !
  !================================================================================================================================
  !

  !>Finalises an interface and deallocates all memory.
  SUBROUTINE INTERFACE_FINALISE(INTERFACE,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("INTERFACE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%COUPLED_MESHES)) DEALLOCATE(INTERFACE%COUPLED_MESHES)
      IF(ASSOCIATED(INTERFACE%MESH_CONNECTIVITY)) THEN
        CALL INTERFACE_MESH_CONNECTIVITY_FINALISE(INTERFACE%MESH_CONNECTIVITY,ERR,ERROR,*999)
      ENDIF
      IF(ASSOCIATED(INTERFACE%POINTS_CONNECTIVITY))THEN
      ENDIF
      IF(ASSOCIATED(INTERFACE%NODES)) CALL NODES_DESTROY(INTERFACE%NODES,ERR,ERROR,*999)
      CALL MESHES_FINALISE(INTERFACE%MESHES,ERR,ERROR,*999)
      CALL FIELDS_FINALISE(INTERFACE%FIELDS,ERR,ERROR,*999)
      CALL INTERFACE_CONDITIONS_FINALISE(INTERFACE%INTERFACE_CONDITIONS,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE)
    ENDIF
    
    CALL EXITS("INTERFACE_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface.
  SUBROUTINE INTERFACE_INITIALISE(INTERFACE,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("INTERFACE_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      CALL FLAG_ERROR("Interface is already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(INTERFACE,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface.",ERR,ERROR,*999)
      INTERFACE%USER_NUMBER=0
      INTERFACE%GLOBAL_NUMBER=0
      INTERFACE%INTERFACE_FINISHED=.FALSE.
      INTERFACE%LABEL=""
      NULLIFY(INTERFACE%INTERFACES)
      NULLIFY(INTERFACE%PARENT_REGION)
      INTERFACE%NUMBER_OF_COUPLED_MESHES=0
      NULLIFY(INTERFACE%COUPLED_MESHES)
      NULLIFY(INTERFACE%MESH_CONNECTIVITY)
      NULLIFY(INTERFACE%POINTS_CONNECTIVITY)
      NULLIFY(INTERFACE%NODES)
      NULLIFY(INTERFACE%MESHES)
      NULLIFY(INTERFACE%GENERATED_MESHES)
      NULLIFY(INTERFACE%FIELDS)
      NULLIFY(INTERFACE%INTERFACE_CONDITIONS)
      NULLIFY(INTERFACE%COORDINATE_SYSTEM)
      NULLIFY(INTERFACE%DATA_POINTS)
      CALL MESHES_INITIALISE(INTERFACE,ERR,ERROR,*999)
      CALL GENERATED_MESHES_INITIALISE(INTERFACE,ERR,ERROR,*999)
      CALL FIELDS_INITIALISE(INTERFACE,ERR,ERROR,*999)
      CALL INTERFACE_CONDITIONS_INITIALISE(INTERFACE,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_INITIALISE")
    RETURN
999 CALL ERRORS("INTERFACE_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the label of an interface for a character label. \see OPENCMISS::CMISSInterfaceLabelGet
  SUBROUTINE INTERFACE_LABEL_GET_C(INTERFACE,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: LABEL !<On return the interface label.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: C_LENGTH,VS_LENGTH

    CALL ENTERS("INTERFACE_LABEL_GET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      C_LENGTH=LEN(LABEL)
      VS_LENGTH=LEN_TRIM(INTERFACE%LABEL)
      IF(C_LENGTH>VS_LENGTH) THEN
        LABEL=CHAR(INTERFACE%LABEL,VS_LENGTH)
      ELSE
        LABEL=CHAR(INTERFACE%LABEL,C_LENGTH)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_LABEL_GET_C")
    RETURN
999 CALL ERRORS("INTERFACE_LABEL_GET_C",ERR,ERROR)
    CALL EXITS("INTERFACE_LABEL_GET_C")
    RETURN 1
    
  END SUBROUTINE INTERFACE_LABEL_GET_C

   !
  !================================================================================================================================
  !

  !>Returns the label of an interface for a varying string label. \see OPENCMISS::CMISSInterfaceLabelGet
  SUBROUTINE INTERFACE_LABEL_GET_VS(INTERFACE,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: LABEL !<On return the interface label.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_LABEL_GET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      !\todo The following line crashes the AIX compiler unless it has a VAR_STR(CHAR()) around it
      LABEL=VAR_STR(CHAR(INTERFACE%LABEL))
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_LABEL_GET_VS")
    RETURN
999 CALL ERRORS("INTERFACE_LABEL_GET_VS",ERR,ERROR)
    CALL EXITS("INTERFACE_LABEL_GET_VS")
    RETURN 1
    
  END SUBROUTINE INTERFACE_LABEL_GET_VS

  !
  !================================================================================================================================
  !

  !>Sets the label of an interface for a character label. \see OPENCMISS::CMISSInterfaceLabelSet
  SUBROUTINE INTERFACE_LABEL_SET_C(INTERFACE,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_LABEL_SET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FLAG_ERROR("Interface has been finished.",ERR,ERROR,*999)
      ELSE
        INTERFACE%LABEL=LABEL
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_LABEL_SET_C")
    RETURN
999 CALL ERRORS("INTERFACE_LABEL_SET_C",ERR,ERROR)
    CALL EXITS("INTERFACE_LABEL_SET_C")
    RETURN 1
  END SUBROUTINE INTERFACE_LABEL_SET_C

  !
  !================================================================================================================================
  !

  !>Sets the label of an interface for a varying string label. \see OPENCMISS::CMISSInterfaceLabelSet
  SUBROUTINE INTERFACE_LABEL_SET_VS(INTERFACE,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_LABEL_SET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FLAG_ERROR("Interface has been finished.",ERR,ERROR,*999)
      ELSE
        INTERFACE%LABEL=LABEL
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_LABEL_SET_VS")
    RETURN
999 CALL ERRORS("INTERFACE_LABEL_SET_VS",ERR,ERROR)
    CALL EXITS("INTERFACE_LABEL_SET_VS")
    RETURN 1
  END SUBROUTINE INTERFACE_LABEL_SET_VS

  !
  !================================================================================================================================
  !

  !>Finalises a meshes connectivity for an interface.
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH(INTERFACE_MESH_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface meshes connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH",ERR,ERROR,*999)

     IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
       IF(INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
         CALL FLAG_ERROR("Interface meshes connectivity has already been finished.",ERR,ERROR,*999)
       ELSE
         !Calculate line or face numbers for coupled mesh elements that are connected to the interface mesh
         CALL INTERFACE_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE(INTERFACE_MESH_CONNECTIVITY,ERR,ERROR,*999)
         INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED=.TRUE.
       ENDIF
     ELSE
       CALL FLAG_ERROR("Interface meshes connectivity is not associated.",ERR,ERROR,*999)
     ENDIF
    
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Initialises a meshes connectivity for an interface.
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CREATE_START(INTERFACE,MESH,INTERFACE_MESH_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the meshes connectivity for
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<On return, a pointer to the created meshes connectivity
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE%MESH_CONNECTIVITY)) THEN
          CALL FLAG_ERROR("The interface already has a meshes connectivity associated.",ERR,ERROR,*999)
        ELSE
          !Initialise the meshes connectivity
          CALL INTERFACE_MESH_CONNECTIVITY_INITIALISE(INTERFACE,MESH,ERR,ERROR,*999)
          !Return the pointer
          INTERFACE_MESH_CONNECTIVITY=>INTERFACE%MESH_CONNECTIVITY
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_CREATE_START")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises a meshes connectivity and deallocates all memory
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_DESTROY(INTERFACE_MESH_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface meshes connectivity to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      CALL INTERFACE_MESH_CONNECTIVITY_FINALISE(INTERFACE_MESH_CONNECTIVITY,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Interface meshes connectivity is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_DESTROY")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_DESTROY",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_DESTROY")
    RETURN 1
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_DESTROY

  !
  !================================================================================================================================
  !

  !>Sets the interface mesh connectivity basis
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_SET_BASIS(INTERFACE_MESH_CONNECTIVITY,BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to interface mesh connectivity to set the element number of elements for.
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_SET_BASIS",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      IF(INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Interface mesh connectivity already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY%BASIS)) THEN
          CALL FLAG_ERROR("Mesh connectivity basis already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(BASIS)) THEN
            INTERFACE_MESH_CONNECTIVITY%BASIS=>BASIS
          ELSE
            CALL FLAG_ERROR("Basis to set mesh connectivity not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mesh connectivity is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_SET_BASIS")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_SET_BASIS",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_SET_BASIS")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_SET_BASIS

  !
  !================================================================================================================================
  !
    
  !>Sets the mapping from an xi position of a coupled mesh element to a node of an interface mesh element
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET(INTERFACE_MESH_CONNECTIVITY,INTERFACE_MESH_ELEMENT_NUMBER, &
    & COUPLED_MESH_INDEX,COUPLED_MESH_ELEMENT_NUMBER,INTERFACE_MESH_LOCAL_NODE_NUMBER,INTERFACE_MESH_COMPONENT_NUMBER,XI, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface mesh connectivity for the interface mesh
    INTEGER(INTG), INTENT(IN) :: INTERFACE_MESH_ELEMENT_NUMBER !<The interface mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_INDEX !<The index of the coupled mesh at the interface to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_ELEMENT_NUMBER !<The coupled mesh element to define the element xi connectivity from
    INTEGER(INTG), INTENT(IN) :: INTERFACE_MESH_LOCAL_NODE_NUMBER !<The interface mesh node to assign the coupled mesh element xi to
    INTEGER(INTG), INTENT(IN) :: INTERFACE_MESH_COMPONENT_NUMBER !<The interface mesh node's component to assign the coupled mesh element xi to
    REAL(DP), INTENT(IN) :: XI(:) !<XI(xi_idx). The xi value for the xi_idx'th xi direction in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    INTEGER(INTG) :: MESH_IDX!<The index of self-contacting mesh, if not a self-contact problem it's the same as the coupled mesh index
    
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      IF(INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Interface mesh connectivity already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
          IF((INTERFACE_MESH_ELEMENT_NUMBER>0).AND. &
            & (INTERFACE_MESH_ELEMENT_NUMBER<=INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS)) THEN
            IF((COUPLED_MESH_INDEX>0).AND.(COUPLED_MESH_INDEX<=INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES)) THEN
              IF((COUPLED_MESH_ELEMENT_NUMBER>0).AND.(COUPLED_MESH_ELEMENT_NUMBER<= &
                & INTERFACE_MESH_CONNECTIVITY%INTERFACE%COUPLED_MESHES(COUPLED_MESH_INDEX)%PTR%NUMBER_OF_ELEMENTS))THEN
                IF((INTERFACE_MESH_COMPONENT_NUMBER>0).AND. &
                  & (INTERFACE_MESH_COMPONENT_NUMBER<=INTERFACE_MESH_CONNECTIVITY%INTERFACE_MESH%NUMBER_OF_COMPONENTS)) THEN
                  IF((INTERFACE_MESH_LOCAL_NODE_NUMBER>0).AND.(INTERFACE_MESH_LOCAL_NODE_NUMBER<= &
                    & INTERFACE_MESH_CONNECTIVITY%BASIS%NUMBER_OF_NODES))THEN
                    IF (INTERFACE_MESH_CONNECTIVITY%INTERFACE%SELF_CONTACT .AND. INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY &
                        & (INTERFACE_MESH_ELEMENT_NUMBER,COUPLED_MESH_INDEX+1)%COUPLED_MESH_ELEMENT_NUMBER .NE. 0) THEN
                      MESH_IDX=COUPLED_MESH_INDEX+1
                    ELSE
                      MESH_IDX=COUPLED_MESH_INDEX
                    ENDIF
                    ELEMENT_CONNECTIVITY=>INTERFACE_MESH_CONNECTIVITY% &
                      & ELEMENT_CONNECTIVITY(INTERFACE_MESH_ELEMENT_NUMBER,MESH_IDX)
                    IF(ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER==COUPLED_MESH_ELEMENT_NUMBER)THEN
                      ELEMENT_CONNECTIVITY%XI(:,INTERFACE_MESH_COMPONENT_NUMBER,INTERFACE_MESH_LOCAL_NODE_NUMBER)=XI(:)
                    ELSE
                      CALL FLAG_ERROR("Coupled mesh element number doesn't match that set to the interface.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Interface local node number is out of range.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Interface component number is out of range.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Coupled mesh element number out of range.",ERR,ERROR,*999)
              END IF
            ELSE
              CALL FLAG_ERROR("Interface coupled mesh index number out of range.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface mesh element number out of range.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface elements connectivity array not allocated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mesh connectivity is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET
  
  !
  !================================================================================================================================
  !
  
  !>Finalises the meshes connectivity and deallocates all memory
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_CONTACT_SET(MESH_CON,INT_ELEM,COUPLED_MESHID,COUPLED_ELEM,LOCAL_CONTACT, &
    & LOCAL_NODE,COMP_NO,XI,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: MESH_CON !<A pointer to interface mesh connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: INT_ELEM !<
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESHID !<The index of the coupled mesh in the interface to set the number of elements for.
    INTEGER(INTG), INTENT(IN) :: COUPLED_ELEM !<
    INTEGER(INTG), INTENT(IN) :: LOCAL_CONTACT !<Local line/face number of the coupled element
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE !<
    INTEGER(INTG), INTENT(IN) :: COMP_NO !<
    REAL(DP), INTENT(IN) :: XI(:) !<XI(xi_idx). The xi value for the xi_idx'th xi direction in the coupled mesh element, can be one dimension less.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    !Local Variables
    INTEGER(INTG) :: I, K !< Dummy index
    INTEGER(INTG) :: xi_dir !<Number of XI directions of the coupled mesh
    INTEGER(INTG) :: MESH_IDX !<The index of self-contacting mesh, if not a self-contact problem it's the same as the coupled mesh index
    
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_CONTACT_SET",ERR,ERROR,*999)

    ! Preliminary error checks to verify user input information
    IF(.NOT.ASSOCIATED(MESH_CON)) CALL FLAG_ERROR("Interface mesh connectivity is not associated.",ERR,ERROR,*999)
    IF(MESH_CON%MESH_CONNECTIVITY_FINISHED) CALL FLAG_ERROR("Interface mesh connectivity already been finished.",ERR,ERROR,*999)
    IF (.NOT.ALLOCATED(MESH_CON%ELEMENT_CONNECTIVITY)) CALL FLAG_ERROR("Interface elements connectivity array not allocated.", &
      & ERR,ERROR,*999)
    IF((INT_ELEM > MESH_CON%NUMBER_OF_INTERFACE_ELEMENTS).OR.(INT_ELEM < 0)) &
      & CALL FLAG_ERROR("Interface mesh element number out of range.",ERR,ERROR,*999)
    IF((COUPLED_MESHID > MESH_CON%NUMBER_OF_COUPLED_MESHES).OR.(COUPLED_MESHID < 0)) &
      & CALL FLAG_ERROR("Interface coupled mesh index number out of range.",ERR,ERROR,*999)
    IF((COUPLED_ELEM>MESH_CON%INTERFACE%COUPLED_MESHES(COUPLED_MESHID)%PTR%NUMBER_OF_ELEMENTS).OR.(COUPLED_ELEM<0))THEN
      CALL FLAG_ERROR("Coupled mesh element number out of range.",ERR,ERROR,*999)
    END IF
    IF((COMP_NO<0).OR.(COMP_NO>MESH_CON%INTERFACE_MESH%NUMBER_OF_COMPONENTS+1)) CALL FLAG_ERROR("Interface component number is &
      & out of range.",ERR,ERROR,*999)
    IF((LOCAL_NODE<0).OR.(LOCAL_NODE> MESH_CON%BASIS%NUMBER_OF_NODES))THEN
      CALL FLAG_ERROR("Interface local node number is out of range.",ERR,ERROR,*999)   
    END IF

    ! Core routine 
    IF(MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,COUPLED_MESHID)%COUPLED_MESH_ELEMENT_NUMBER/=COUPLED_ELEM)THEN
      CALL FLAG_ERROR("Coupled mesh element number doesn't match that set to the interface.",ERR,ERROR,*999)
    ELSE
      IF (MESH_CON%INTERFACE%SELF_CONTACT .AND. MESH_CON%ELEMENT_CONNECTIVITY &
          & (INT_ELEM,COUPLED_MESHID+1)%COUPLED_MESH_ELEMENT_NUMBER .NE. 0) THEN
        MESH_IDX=COUPLED_MESHID+1
      ELSE
        MESH_IDX=COUPLED_MESHID
      ENDIF
      IF (SIZE(XI,1)==MESH_CON%BASIS%NUMBER_OF_XI) THEN
        SELECT CASE(MESH_CON%BASIS%NUMBER_OF_XI)
        CASE(1)
          MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%CONNECTED_LINE=LOCAL_CONTACT
          SELECT CASE(LOCAL_CONTACT)
          CASE(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=0.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=-2
          CASE(2)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=1.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=2
          CASE(3)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=0.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=-1
          CASE(4)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=1.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=1
          CASE DEFAULT 
            LOCAL_ERROR="The local contact number for coupled mesh "// &
                & TRIM(NUMBER_TO_VSTRING(MESH_IDX, &
                & "*",ERR,ERROR))//" is invalid for element number "//TRIM(NUMBER_TO_VSTRING(COUPLED_ELEM,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(2)
          MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%CONNECTED_FACE=LOCAL_CONTACT
          SELECT CASE(LOCAL_CONTACT)
          CASE(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=1.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(3,COMP_NO,LOCAL_NODE)=XI(2)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=1
          CASE(2)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=0.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(3,COMP_NO,LOCAL_NODE)=XI(2)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=-1
          CASE(3)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=1.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(3,COMP_NO,LOCAL_NODE)=XI(2)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=2
          CASE(4)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=0.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(3,COMP_NO,LOCAL_NODE)=XI(2)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=-2
          CASE(5)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(3,COMP_NO,LOCAL_NODE)=1.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=XI(2)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=3
          CASE(6)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(3,COMP_NO,LOCAL_NODE)=0.0_DP
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(1,COMP_NO,LOCAL_NODE)=XI(1)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(2,COMP_NO,LOCAL_NODE)=XI(2)
            MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%COUPLED_MESH_CONTACT_XI_NORMAL=-3
          CASE DEFAULT 
            LOCAL_ERROR="The local contact number for coupled mesh "// &
                & TRIM(NUMBER_TO_VSTRING(COUPLED_MESHID, &
                & "*",ERR,ERROR))//" is invalid for element number "//TRIM(NUMBER_TO_VSTRING(COUPLED_ELEM,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="interface dimension should be < 3"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        MESH_CON%ELEMENT_CONNECTIVITY(INT_ELEM,MESH_IDX)%XI(:,COMP_NO,LOCAL_NODE)=XI(:)
      END IF
    END IF
    
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_CONTACT_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_CONTACT_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_CONTACT_SET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_CONTACT_SET

  !
  !================================================================================================================================
  !
  
  !>Sets the connectivity between an element in a coupled mesh to an element in the interface mesh
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET(INTERFACE_MESH_CONNECTIVITY,INTERFACE_MESH_ELEMENT_NUMBER, &
      & COUPLED_MESH_INDEX,COUPLED_MESH_ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface mesh connectivity for the interface mesh
    INTEGER(INTG), INTENT(IN) :: INTERFACE_MESH_ELEMENT_NUMBER !<The interface mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_INDEX !<The index of the coupled mesh at the interface to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_ELEMENT_NUMBER !<The coupled mesh element to be connected to the interface
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NumberOfInterfaceElementNodes,NumberOfCoupledMeshXiDirections
    INTEGER(INTG) :: MESH_IDX!<The index of self-contacting mesh, if not a self-contact problem it's the same as the coupled mesh index
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      IF(INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Interface mesh connectivity has already been finished.",ERR,ERROR,*999)
      ELSE
        IF((INTERFACE_MESH_ELEMENT_NUMBER>0).AND.(INTERFACE_MESH_ELEMENT_NUMBER<= &
          & INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS)) THEN
          IF((COUPLED_MESH_INDEX>0).AND.(COUPLED_MESH_INDEX<=INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES)) THEN
            IF (ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
              IF (INTERFACE_MESH_CONNECTIVITY%INTERFACE%SELF_CONTACT .AND. INTERFACE_MESH_CONNECTIVITY% &
                  & ELEMENT_CONNECTIVITY(INTERFACE_MESH_ELEMENT_NUMBER,COUPLED_MESH_INDEX)%COUPLED_MESH_ELEMENT_NUMBER .NE. 0) THEN
                MESH_IDX=COUPLED_MESH_INDEX+1
              ELSE
                MESH_IDX=COUPLED_MESH_INDEX
              ENDIF
              ELEMENT_CONNECTIVITY=>INTERFACE_MESH_CONNECTIVITY% &
                & ELEMENT_CONNECTIVITY(INTERFACE_MESH_ELEMENT_NUMBER,MESH_IDX)
              ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER=COUPLED_MESH_ELEMENT_NUMBER
              !\todo NumberOfCoupledMeshXiDirections currently set to the number of interface mesh xi directions + 1. Restructure ELEMENT_CONNECTIVITY type see below
              NumberOfCoupledMeshXiDirections=INTERFACE_MESH_CONNECTIVITY%INTERFACE_MESH%NUMBER_OF_DIMENSIONS+1
              NumberOfInterfaceElementNodes=INTERFACE_MESH_CONNECTIVITY%BASIS%NUMBER_OF_NODES
              IF(ALLOCATED(ELEMENT_CONNECTIVITY%XI)) THEN
                LOCAL_ERROR="Interface mesh element connectivity already allocated for coupled mesh element " &
                  & //TRIM(NUMBER_TO_VSTRING(COUPLED_MESH_ELEMENT_NUMBER,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                !\todo Update mesh component index to look at the number of mesh components in each element. 
                !\todo Currently this defaults to the first mesh component ie %XI(NumberOfInterfaceMeshXi,1,NumberOfInterfaceElementNodes)). 
                !\todo The interface mesh types will also need to be restructured.
                !eg. ELEMENT_CONNECTIVITY%COUPLED_MESH(CoupledMeshIdx)%MESH_COMPONENT(MeshComponentIdx)%XI(NumberOfCoupledMeshXiDirections,NumberOfInterfaceElementNodes) and adding appropriate initialize/finialize routines
                ALLOCATE(ELEMENT_CONNECTIVITY%XI(NumberOfCoupledMeshXiDirections,1,NumberOfInterfaceElementNodes),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface element connectivity.",ERR,ERROR,*999)
                ELEMENT_CONNECTIVITY%XI=0.0_DP
                ELEMENT_CONNECTIVITY%CONNECTED_LINE=0
                ELEMENT_CONNECTIVITY%CONNECTED_FACE=0
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface elements connectivity array not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface coupled mesh index number out of range.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mesh element number out of range.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mesh connectivity is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET

  !
  !================================================================================================================================
  !
  
  !>Calculate line or face numbers for coupled mesh elements that are connected to the interface mesh
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE(INTERFACE_MESH_CONNECTIVITY,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to interface mesh connectivity to calculate line or face numbers for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: XiDifference
    INTEGER(INTG) :: CoupledMeshDomainIdx,InterfaceElementIdx,CoupledMeshXiIdx,NumberOfInterfaceElementNodes, &
      & NumberOfInterfaceMeshXi,CoupledMeshIdx
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE",ERR,ERROR,*999)
    
    DO CoupledMeshDomainIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
      DO InterfaceElementIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS
        IF(INTERFACE_MESH_CONNECTIVITY%INTERFACE%SELF_CONTACT) THEN
          CoupledMeshIdx=1 !Only 1 mesh in self_contact case
        ELSE
          CoupledMeshIdx=CoupledMeshDomainIdx
        ENDIF
        ELEMENT_CONNECTIVITY=>INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)
        NumberOfInterfaceElementNodes=INTERFACE_MESH_CONNECTIVITY%BASIS%NUMBER_OF_NODES
        NumberOfInterfaceMeshXi=INTERFACE_MESH_CONNECTIVITY%INTERFACE_MESH%NUMBER_OF_DIMENSIONS
        SELECTCASE(NumberOfInterfaceMeshXi)
        CASE(1) !Lines
          DO CoupledMeshXiIdx=1,INTERFACE_MESH_CONNECTIVITY%INTERFACE%COUPLED_MESHES(CoupledMeshIdx)%PTR%NUMBER_OF_DIMENSIONS
            ! Calculate difference between first node and last node of an element
            XiDifference=ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,1)- &
              & ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,NumberOfInterfaceElementNodes)
            IF(XiDifference==0.0_DP) THEN
              IF(ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,NumberOfInterfaceElementNodes)==0.0_DP) THEN
                ELEMENT_CONNECTIVITY%CONNECTED_LINE=3-(CoupledMeshXiIdx-1)*2
              ELSE
                ELEMENT_CONNECTIVITY%CONNECTED_LINE=4-(CoupledMeshXiIdx-1)*2
              ENDIF
            ENDIF
          ENDDO
          SELECT CASE(ELEMENT_CONNECTIVITY%CONNECTED_LINE)
          CASE(1)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=-2
          CASE(2)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=2
          CASE(3)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=-1
          CASE(4)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=1
          CASE DEFAULT
            LOCAL_ERROR="Local line number should be < 4"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(2) !Faces
          DO CoupledMeshXiIdx=1,INTERFACE_MESH_CONNECTIVITY%INTERFACE%COUPLED_MESHES(CoupledMeshIdx)%PTR%NUMBER_OF_DIMENSIONS
            XiDifference=ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,1)- &
              & ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,NumberOfInterfaceElementNodes)
            IF(XiDifference==0.0_DP) THEN
              IF(ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,NumberOfInterfaceElementNodes)==0.0_DP) THEN
                ELEMENT_CONNECTIVITY%CONNECTED_FACE=(CoupledMeshXiIdx-1)*2+2
              ELSE
                ELEMENT_CONNECTIVITY%CONNECTED_FACE=(CoupledMeshXiIdx-1)*2+1
              ENDIF
            ENDIF
          ENDDO
          SELECT CASE(ELEMENT_CONNECTIVITY%CONNECTED_FACE)
          CASE(1)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=1
          CASE(2)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=-1
          CASE(3)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=2
          CASE(4)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=-2
          CASE(5)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=3
          CASE(6)
            ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL=-3  
          CASE DEFAULT
            LOCAL_ERROR="Local face number should be < 6"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT  
        CASE DEFAULT 
          LOCAL_ERROR="Number of interface mesh dimension of "//TRIM(NUMBER_TO_VSTRING(NumberOfInterfaceMeshXi,"*",ERR,ERROR))// &
            & "is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDSELECT
      ENDDO
    ENDDO
    
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE")
    RETURN
    
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE")
    RETURN 1
  
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the meshes connectivity and deallocates all memory
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_FINALISE(INTERFACE_MESH_CONNECTIVITY,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE) :: INTERFACE_MESH_CONNECTIVITY !<The interface mesh connectivity to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_FINALISE",ERR,ERROR,*999)

    CALL INTERFACE_MESH_CONNECTIVITY_ELEMENT_FINALISE(INTERFACE_MESH_CONNECTIVITY,ERR,ERROR,*999)
    NULLIFY(INTERFACE_MESH_CONNECTIVITY%INTERFACE)
    NULLIFY(INTERFACE_MESH_CONNECTIVITY%INTERFACE_MESH)
    NULLIFY(INTERFACE_MESH_CONNECTIVITY%BASIS)
    INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS=0
    INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES=0
       
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_FINALISE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface mesh connectivity.
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_INITIALISE(INTERFACE,MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the mesh connectivity for
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%MESH_CONNECTIVITY)) THEN
        CALL FLAG_ERROR("Interface mesh connectivity is already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE%MESH_CONNECTIVITY,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface mesh connectivity.",ERR,ERROR,*999)
        INTERFACE%MESH_CONNECTIVITY%INTERFACE=>INTERFACE
        INTERFACE%MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED=.FALSE.
        INTERFACE%MESH_CONNECTIVITY%INTERFACE_MESH=>MESH
        NULLIFY(INTERFACE%MESH_CONNECTIVITY%BASIS)
        CALL INTERFACE_MESH_CONNECTIVITY_ELEMENT_INITIALISE(INTERFACE,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_INITIALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the meshes connectivity and deallocates all memory
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_FINALISE(INTERFACE_POINTS_CONNECTIVITY,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE) :: INTERFACE_POINTS_CONNECTIVITY !<The interface points connectivity to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_FINALISE",ERR,ERROR,*999)

    CALL INTERFACE_POINTS_CONNECTIVITY_POINTS_FINALISE(INTERFACE_POINTS_CONNECTIVITY,ERR,ERROR,*999)
    NULLIFY(INTERFACE_POINTS_CONNECTIVITY%INTERFACE)
    NULLIFY(INTERFACE_POINTS_CONNECTIVITY%INTERFACE_MESH)
    INTERFACE_POINTS_CONNECTIVITY%NUMBER_OF_ELEMENTS=0
    INTERFACE_POINTS_CONNECTIVITY%NUMBER_OF_DATA_POINTS=0
    INTERFACE_POINTS_CONNECTIVITY%NUMBER_INT_DOM=0
       
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_FINALISE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface data points connectivity.
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_INITIALISE(INTERFACE,MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the mesh connectivity for
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("INTERFACE_POINT_CONNECTIVITY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%POINTS_CONNECTIVITY)) THEN
        CALL FLAG_ERROR("Interface points connectivity is already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE%POINTS_CONNECTIVITY,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface points connectivity.",ERR,ERROR,*999)
        INTERFACE%POINTS_CONNECTIVITY%INTERFACE=>INTERFACE
        INTERFACE%POINTS_CONNECTIVITY%POINTS_CONNECTIVITY_FINISHED=.FALSE.
        INTERFACE%POINTS_CONNECTIVITY%INTERFACE_MESH=>MESH
        CALL INTERFACE_POINTS_CONNECTIVITY_POINTS_INITIALISE(INTERFACE,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_INITIALISE")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises data points connectivity for an interface.
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_CREATE_START(INTERFACE,MESH,INTERFACE_POINTS_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the meshes connectivity for
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: INTERFACE_POINTS_CONNECTIVITY !<On return, a pointer to the created points connectivity
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE%POINTS_CONNECTIVITY)) THEN
          CALL FLAG_ERROR("The interface already has a points connectivity associated.",ERR,ERROR,*999)
        ELSE
          !Initialise the poins connectivity
          CALL INTERFACE_POINTS_CONNECTIVITY_INITIALISE(INTERFACE,MESH,ERR,ERROR,*999)
          !Return the pointer
          INTERFACE_POINTS_CONNECTIVITY=>INTERFACE%POINTS_CONNECTIVITY
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_CREATE_START")
    RETURN 1
    
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_CREATE_START
  
  !
  !================================================================================================================================
  !
  
  !>Finalises the meshes connectivity and deallocates all memory
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_ELEMENT_NUMBER_SET(POINTS_CON,POINTS_IDX,COUPLED_MESHID,COUPLED_ELEM,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: POINTS_CON !<A pointer to interface points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: POINTS_IDX !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESHID !<The index of the coupled mesh in the interface to set the number of elements for.
    INTEGER(INTG), INTENT(IN) :: COUPLED_ELEM !<The coupled mesh element number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: xi_dir !<Number of XI directions of the coupled mesh
    
    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_ELEMENT_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(POINTS_CON)) THEN
      IF(POINTS_CON%POINTS_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Interface points connectivity has already been finished.",ERR,ERROR,*999)
      ELSE
        IF ((POINTS_IDX > POINTS_CON%NUMBER_OF_DATA_POINTS).OR.(POINTS_IDX < 0)) THEN
          CALL FLAG_ERROR("Interface data points index number out of range.",ERR,ERROR,*999)
        ELSE
          IF ((COUPLED_MESHID > POINTS_CON%NUMBER_INT_DOM).OR.(COUPLED_MESHID < 0)) THEN
            CALL FLAG_ERROR("Interface coupled mesh index number out of range.",ERR,ERROR,*999)
          ELSE
            IF (ALLOCATED(POINTS_CON%POINTS_CONNECTIVITY)) THEN
              POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_ELEMENT_NUMBER=COUPLED_ELEM
              xi_dir=POINTS_CON%INTERFACE_MESH%NUMBER_OF_DIMENSIONS+1
              IF(.NOT.ALLOCATED(POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI)) THEN
                !\todo Update mesh component index to look at the number of mesh components in each element. Currently this defaults to the first mesh component ie %XI(xi_dir,1)). The interface mesh types will also need to be restructured.
                ALLOCATE(POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(xi_dir,1))
              ENDIF
              POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI=0.5_DP
            ELSE
              CALL FLAG_ERROR("Interface points connectivity array not allocated.",ERR,ERROR,*999)
            END IF
          END IF
        END IF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface points connectivity is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_ELEMENT_NUMBER_SET")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_ELEMENT_NUMBER_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_ELEMENT_NUMBER_SET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_ELEMENT_NUMBER_SET
  
  !
  !================================================================================================================================
  !
    
  !>Finalises the meshes connectivity and deallocates all memory
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_POINT_XI_SET(POINTS_CON,POINTS_IDX,COUPLED_MESHID,COUPLED_ELEM, &
    & COMP_NO,XI,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: POINTS_CON !<A pointer to interface points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: POINTS_IDX !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESHID !<The index of the coupled mesh in the interface to set the number of elements for.
    INTEGER(INTG), INTENT(IN) :: COUPLED_ELEM !<The coupled mesh element number
    INTEGER(INTG), INTENT(IN) :: COMP_NO !<Mesh component number
    REAL(DP), INTENT(IN) :: XI(:) !<XI(xi_idx). The xi value for the xi_idx'th xi direction in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_POINT_XI_SET",ERR,ERROR,*999)

    ! Preliminary error checks to verify user input information
    IF(ASSOCIATED(POINTS_CON)) THEN
      IF(POINTS_CON%POINTS_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Interface mesh connectivity already been finished.",ERR,ERROR,*999)
      ELSE
        IF (ALLOCATED(POINTS_CON%POINTS_CONNECTIVITY)) THEN
          IF((POINTS_IDX > 0).OR.(POINTS_IDX < POINTS_CON%NUMBER_OF_DATA_POINTS)) THEN
            IF((COUPLED_MESHID > 0).OR.(COUPLED_MESHID < POINTS_CON%NUMBER_INT_DOM)) THEN
              IF((COUPLED_ELEM>0).OR.(COUPLED_ELEM<POINTS_CON%INTERFACE%COUPLED_MESHES(COUPLED_MESHID)%PTR%NUMBER_OF_ELEMENTS))THEN
                IF((COMP_NO>0).OR.(COMP_NO<POINTS_CON%INTERFACE_MESH%NUMBER_OF_COMPONENTS+1)) THEN         
                  !Core routine 
                  IF(POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_ELEMENT_NUMBER/=COUPLED_ELEM)THEN
                    CALL FLAG_ERROR("Coupled mesh element number doesn't match that set to the interface.",ERR,ERROR,*999)
                  ELSE
                    POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(:,COMP_NO)=XI(:) 
                  END IF
                ELSE
                  CALL FLAG_ERROR("Interface component number is out of range.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Coupled mesh element number out of range.",ERR,ERROR,*999)
              END IF
            ELSE
              CALL FLAG_ERROR("Interface coupled mesh index number out of range.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface data point index out of range.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface elements connectivity array not allocated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface points connectivity is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_POINT_XI_SET")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_POINT_XI_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_POINT_XI_SET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_POINT_XI_SET
  
  !
  !================================================================================================================================
  !

  !>Initialises a meshes connectivity for an interface.
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_CREATE_FINISH(INTERFACE_POINTS_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: INTERFACE_POINTS_CONNECTIVITY !<A pointer to the interface poiints connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx
  
    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_CREATE_FINISH",ERR,ERROR,*999)

     IF(ASSOCIATED(INTERFACE_POINTS_CONNECTIVITY)) THEN
       IF(INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY_FINISHED) THEN
         CALL FLAG_ERROR("Interface points connectivity has already been finished.",ERR,ERROR,*999)
       ELSE
         IF(INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(1,1)%COUPLED_MESH_CONTACT_NUMBER==0) THEN
           CALL INTERFACE_POINTS_CONNECTIVITY_CONTACT_NUMBER_CALCULATE(INTERFACE_POINTS_CONNECTIVITY, &
             & ERR,ERROR,*999)!also calculate normal xi direction
         ENDIF
         INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY_FINISHED=.TRUE.
         INTERFACE_POINTS_CONNECTIVITY%NUMBER_OF_ELEMENTS=INTERFACE_POINTS_CONNECTIVITY%INTERFACE%MESHES%MESHES(1)% &
           & PTR%NUMBER_OF_ELEMENTS
         ALLOCATE(INTERFACE_POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS(INTERFACE_POINTS_CONNECTIVITY%NUMBER_OF_ELEMENTS, &
           & INTERFACE_POINTS_CONNECTIVITY%NUMBER_INT_DOM),STAT=ERR)
       ENDIF
     ELSE
       CALL FLAG_ERROR("Interface points connectivity is not associated.",ERR,ERROR,*999)
     ENDIF
    
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_CREATE_FINISH
  
  !
  !================================================================================================================================
  !
  
  !>Calculate line/face numbers for coupled mesh elements that are connected to the interface mesh
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_CONTACT_NUMBER_CALCULATE(POINTS_CON,ERR,ERROR,*)
    !Argument variables
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: POINTS_CON !<A pointer to interface points connectivity to calculate line numbers for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Local variables
    REAL(DP) :: XI_DIFF
    INTEGER(INTG) :: COUPLED_MESHID
    INTEGER(INTG) :: point_idx,reference_point_idx,mod_reference_point_idx
    INTEGER(INTG) :: xi_dir
    INTEGER(INTG) :: LOCAL_CONTACT
    LOGICAL :: DIFFERENT_ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_CONTACT_NUMBER_CALCULATE",ERR,ERROR,*999)
    
    DO COUPLED_MESHID=1,POINTS_CON%NUMBER_INT_DOM
      DO point_idx=1,POINTS_CON%NUMBER_OF_DATA_POINTS
        DIFFERENT_ELEMENT=.TRUE.
        reference_point_idx=point_idx+1
        SELECT CASE(SIZE(POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%XI,1))
        CASE(2)      
          DO xi_dir=1,SIZE(POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%XI,1)    
            DO WHILE (DIFFERENT_ELEMENT) 
              mod_reference_point_idx=MOD(reference_point_idx-1,POINTS_CON%NUMBER_OF_DATA_POINTS)+1
              IF ((POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_ELEMENT_NUMBER == &
                  & POINTS_CON%POINTS_CONNECTIVITY(mod_reference_point_idx,COUPLED_MESHID)%COUPLED_MESH_ELEMENT_NUMBER) &
                  &  .AND. (mod_reference_point_idx/=point_idx)) THEN
                DIFFERENT_ELEMENT=.FALSE.
              ELSE
                IF (MOD(reference_point_idx,POINTS_CON%NUMBER_OF_DATA_POINTS)==point_idx) THEN
                  LOCAL_ERROR="Only one data point in this element"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  reference_point_idx=reference_point_idx+1
                ENDIF
              ENDIF
            ENDDO !(DIFFERENT_ELEMENT) 
            XI_DIFF=POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%XI(xi_dir,1)- &
              & POINTS_CON%POINTS_CONNECTIVITY(mod_reference_point_idx,COUPLED_MESHID)%XI(xi_dir,1)
            IF(XI_DIFF==0.0_DP) THEN !Data points have same xi location (current xi direction)
              !Check if the absolute value is zero 
              IF(POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%XI(xi_dir,1) &
                & ==0.0_DP) THEN
                LOCAL_CONTACT=3-(xi_dir-1)*2 !Line 1 or 3
              ELSE
                LOCAL_CONTACT=4-(xi_dir-1)*2 !Line 2 or 4
              END IF
            END IF
          ENDDO !xi_dir
          
          POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_NUMBER=LOCAL_CONTACT
          SELECT CASE(LOCAL_CONTACT)
          CASE(1)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-2
          CASE(2)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=2
          CASE(3)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-1
          CASE(4)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=1
          CASE DEFAULT
            LOCAL_ERROR="Local line number should be < 4"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(3)
          DO xi_dir=1,SIZE(POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%XI,1)-1
            IF(point_idx==POINTS_CON%NUMBER_OF_DATA_POINTS) THEN
              XI_DIFF=POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%XI(xi_dir,1)- &
                & POINTS_CON%POINTS_CONNECTIVITY(point_idx-1,COUPLED_MESHID)%XI(xi_dir,1)
            ELSE
              XI_DIFF=POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%XI(xi_dir,1)- &
                & POINTS_CON%POINTS_CONNECTIVITY(point_idx+1,COUPLED_MESHID)%XI(xi_dir,1)
            ENDIF
            IF(XI_DIFF==0.0_DP) THEN
              IF(POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%XI(xi_dir,1) &
                & ==0.0_DP) THEN
                LOCAL_CONTACT=(xi_dir-1)*2+2
              ELSE
                LOCAL_CONTACT=(xi_dir-1)*2+1
              END IF
            END IF
          ENDDO
          POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_NUMBER=LOCAL_CONTACT
          SELECT CASE(LOCAL_CONTACT)
          CASE(1)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=1
          CASE(2)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-1
          CASE(3)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=2
          CASE(4)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-2
          CASE(5)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=3
          CASE(6)
            POINTS_CON%POINTS_CONNECTIVITY(point_idx,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-3  
          CASE DEFAULT
            LOCAL_ERROR="Local face number should be < 6"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT  
        CASE DEFAULT 
          LOCAL_ERROR="Mesh dimension should be < 3 & > 1"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDDO
    ENDDO
    
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_CONTACT_NUMBER_CALCULATE")
    RETURN
    
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_CONTACT_NUMBER_CALCULATE",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_CONTACT_NUMBER_CALCULATE")
    RETURN 1
  
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_CONTACT_NUMBER_CALCULATE
  
  !
  !================================================================================================================================
  !
  
  !>Finalises the meshes connectivity and deallocates all memory
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_POINT_XI_CONTACT_SET(POINTS_CON,POINTS_IDX,COUPLED_MESHID,COUPLED_ELEM, &
      & LOCAL_CONTACT,COMP_NO,XI,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: POINTS_CON !<A pointer to interface point connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: POINTS_IDX !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESHID !<The index of the coupled mesh in the interface to set the number of elements for.
    INTEGER(INTG), INTENT(IN) :: COUPLED_ELEM !<The coupled mesh element number
    INTEGER(INTG), INTENT(IN) :: LOCAL_CONTACT !<Local line/face number of the coupled element
    INTEGER(INTG), INTENT(IN) :: COMP_NO !<Mesh component number
    REAL(DP), INTENT(IN) :: XI(:) !<XI(xi_idx). The xi value for the xi_idx'th xi direction in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Local Variables
    INTEGER(INTG) :: xi_dir !<Number of XI directions of the coupled mesh
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_POINT_XI_CONTACT_SET",ERR,ERROR,*999)

    ! Preliminary error checks to verify user input information
    IF(ASSOCIATED(POINTS_CON)) THEN
      IF(POINTS_CON%POINTS_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Interface points connectivity already been finished.",ERR,ERROR,*999)
      ELSE
        IF (ALLOCATED(POINTS_CON%POINTS_CONNECTIVITY)) THEN
          IF((POINTS_IDX > POINTS_CON%NUMBER_OF_DATA_POINTS).OR.(POINTS_IDX < 0)) THEN
            CALL FLAG_ERROR("Interface data points index number out of range.",ERR,ERROR,*999)
          ELSE
            IF((COUPLED_MESHID > POINTS_CON%NUMBER_INT_DOM).OR.(COUPLED_MESHID < 0)) THEN
              CALL FLAG_ERROR("Interface coupled mesh index number out of range.",ERR,ERROR,*999)
            ELSE
              IF((COUPLED_ELEM>POINTS_CON%INTERFACE%COUPLED_MESHES &
                  & (COUPLED_MESHID)%PTR%NUMBER_OF_ELEMENTS).OR.(COUPLED_ELEM<0)) THEN 
                CALL FLAG_ERROR("Coupled mesh element number out of range.",ERR,ERROR,*999)
              ELSE
                IF((COMP_NO<0).OR.(COMP_NO>POINTS_CON%INTERFACE_MESH%NUMBER_OF_COMPONENTS+1)) THEN 
                  CALL FLAG_ERROR("Interface component number is out of range.",ERR,ERROR,*999)
                ELSE
                  ! Core routine 
                  IF(POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_ELEMENT_NUMBER/=COUPLED_ELEM)THEN
                    CALL FLAG_ERROR("Coupled mesh element number doesn't match that set to the interface.",ERR,ERROR,*999)
                  ELSE
                    POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_NUMBER=LOCAL_CONTACT
                    IF (SIZE(XI,1)/=SIZE(POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI,1)) THEN
                      SELECT CASE(SIZE(XI,1))
                      CASE(1)
                        SELECT CASE(LOCAL_CONTACT)
                        CASE(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=0.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-2
                        CASE(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=1.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=2
                        CASE(3)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=0.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-1
                        CASE(4)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=1.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=1
                        CASE DEFAULT 
                          LOCAL_ERROR="The local contact number for coupled mesh "//TRIM(NUMBER_TO_VSTRING(COUPLED_MESHID, &
                            & "*",ERR,ERROR))//" is invalid for element number "//TRIM(NUMBER_TO_VSTRING(COUPLED_ELEM, &
                            & "*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE(2)
                        SELECT CASE(LOCAL_CONTACT)
                        CASE(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=1.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(3,COMP_NO)=XI(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=1
                        CASE(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=0.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(3,COMP_NO)=XI(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-1
                        CASE(3)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=1.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(3,COMP_NO)=XI(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=2
                        CASE(4)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=0.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(3,COMP_NO)=XI(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-2
                        CASE(5)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(3,COMP_NO)=1.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=XI(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=3
                        CASE(6)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(3,COMP_NO)=0.0_DP
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(1,COMP_NO)=XI(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(2,COMP_NO)=XI(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-3
                        CASE DEFAULT 
                          LOCAL_ERROR="The local contact number for coupled mesh "// &
                              & TRIM(NUMBER_TO_VSTRING(COUPLED_MESHID,"*",ERR,ERROR))//" is invalid for element number "// &
                              & TRIM(NUMBER_TO_VSTRING(COUPLED_ELEM,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE DEFAULT
                        LOCAL_ERROR="interface dimension should be < 3"
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%XI(:,COMP_NO)=XI(:)
                      SELECT CASE(SIZE(XI,1))
                      CASE(2)
                        SELECT CASE(LOCAL_CONTACT)
                        CASE(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-2
                        CASE(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=2
                        CASE(3)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-1
                        CASE(4)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=1
                        CASE DEFAULT 
                          LOCAL_ERROR="The local contact number for coupled mesh "//TRIM(NUMBER_TO_VSTRING(COUPLED_MESHID, &
                            & "*",ERR,ERROR))//" is invalid for element number "//TRIM(NUMBER_TO_VSTRING(COUPLED_ELEM, &
                            & "*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE(3)
                        SELECT CASE(LOCAL_CONTACT)
                        CASE(1)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=1
                        CASE(2)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-1
                        CASE(3)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=2
                        CASE(4)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-2
                        CASE(5)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=3
                        CASE(6)
                          POINTS_CON%POINTS_CONNECTIVITY(POINTS_IDX,COUPLED_MESHID)%COUPLED_MESH_CONTACT_XI_NORMAL=-3
                        CASE DEFAULT 
                          LOCAL_ERROR="The local contact number for coupled mesh "// &
                              & TRIM(NUMBER_TO_VSTRING(COUPLED_MESHID,"*",ERR,ERROR))//" is invalid for element number "// &
                              & TRIM(NUMBER_TO_VSTRING(COUPLED_ELEM,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999) 
                        END SELECT
                      CASE DEFAULT 
                        LOCAL_ERROR="interface dimension should be < 3"
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ENDIF
                  ENDIF !end of core code   
                ENDIF    
              ENDIF    
            ENDIF
          ENDIF    
        ELSE
          CALL FLAG_ERROR("Interface points connectivity array not allocated.",ERR,ERROR,*999)
        ENDIF  
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface points connectivity is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_POINT_XI_CONTACT_SET")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_POINT_XI_CONTACT_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_POINT_XI_CONTACT_SET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_POINT_XI_CONTACT_SET

  !
  !================================================================================================================================
  !
  
    !>Set interface points connectivity
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_PROJECTION_RESULTS_SET(POINTS_CON,DATA_PROJECTION,COUPLED_MESH_IDX,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: POINTS_CON !<A pointer to interface points connectivity to set the element number of elements for.
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection that contains projection results
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_IDX
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Local Variables
    TYPE(DATA_PROJECTION_RESULT_TYPE), POINTER :: DATA_PROJECTION_RESULT
    INTEGER(INTG) :: COMP_NO=1 !\TODO:mesh component number is now hard-coded to be one, need to be updated
    INTEGER(INTG) :: data_projection_idx,data_point_idx
    LOGICAL:: POINTS_CON_FINISHED
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_PROJECTION_RESULTS_SET",ERR,ERROR,*999)

    ! Preliminary error checks to verify user input information
    IF(ASSOCIATED(POINTS_CON)) THEN
      IF (ALLOCATED(POINTS_CON%POINTS_CONNECTIVITY)) THEN
        IF(ASSOCIATED(DATA_PROJECTION)) THEN
          POINTS_CON%POINTS_CONNECTIVITY_FINISHED=.FALSE.
          IF(DATA_PROJECTION%DATA_PROJECTION_PROJECTED) THEN
            DO data_point_idx=1,POINTS_CON%NUMBER_OF_DATA_POINTS
              DATA_PROJECTION_RESULT=>DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)
              CALL INTERFACE_POINTS_CONNECTIVITY_ELEMENT_NUMBER_SET(POINTS_CON,data_point_idx,COUPLED_MESH_IDX, &
                & DATA_PROJECTION_RESULT%ELEMENT_NUMBER,ERR,ERROR,*999)
              IF(DATA_PROJECTION_RESULT%ELEMENT_LINE_NUMBER/=0) THEN
                CALL INTERFACE_POINTS_CONNECTIVITY_POINT_XI_CONTACT_SET(POINTS_CON,data_point_idx,COUPLED_MESH_IDX, &
                  & DATA_PROJECTION_RESULT%ELEMENT_NUMBER,DATA_PROJECTION_RESULT%ELEMENT_LINE_NUMBER,COMP_NO, &
                  & DATA_PROJECTION_RESULT%XI,ERR,ERROR,*999)
              ELSEIF(DATA_PROJECTION_RESULT%ELEMENT_FACE_NUMBER/=0) THEN
                CALL INTERFACE_POINTS_CONNECTIVITY_POINT_XI_CONTACT_SET(POINTS_CON,data_point_idx,COUPLED_MESH_IDX, &
                  & DATA_PROJECTION_RESULT%ELEMENT_NUMBER,DATA_PROJECTION_RESULT%ELEMENT_FACE_NUMBER,COMP_NO, &
                  & DATA_PROJECTION_RESULT%XI,ERR,ERROR,*999)
              ENDIF
            ENDDO !data_point_idx      
          ELSE
            LOCAL_ERROR="Data Projection  "//TRIM(NUMBER_TO_VSTRING(data_projection_idx,"*",ERR,ERROR))//" is not projected."        
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF                           
        ELSE
          CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
        ENDIF        
      ELSE
        CALL FLAG_ERROR("Interface points connectivity array not allocated.",ERR,ERROR,*999)
      ENDIF  
    ELSE
      CALL FLAG_ERROR("Interface points connectivity is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_PROJECTION_RESULTS_SET")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_PROJECTION_RESULTS_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_PROJECTION_RESULTS_SET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_PROJECTION_RESULTS_SET

  !
  !================================================================================================================================
  !

  !>Finds and returns in INTERFACE a pointer to the interface identified by USER_NUMBER in the given PARENT_REGION. If no interface with that USER_NUMBER exists INTERFACE is left nullified.
  SUBROUTINE INTERFACE_USER_NUMBER_FIND(USER_NUMBER,PARENT_REGION,INTERFACE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find.
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<The parent region to find the interface in    
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<On return a pointer to the interface with the given user number. If no interface with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(PARENT_REGION)) THEN
      IF(ASSOCIATED(INTERFACE)) THEN
        CALL FLAG_ERROR("Interface is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(INTERFACE)
        IF(ASSOCIATED(PARENT_REGION%INTERFACES)) THEN
          interface_idx=1
          DO WHILE(interface_idx<=PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES.AND..NOT.ASSOCIATED(INTERFACE))
            IF(PARENT_REGION%INTERFACES%INTERFACES(interface_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
              INTERFACE=>PARENT_REGION%INTERFACES%INTERFACES(interface_idx)%PTR
            ELSE
              interface_idx=interface_idx+1
            ENDIF
          ENDDO
        ELSE
          LOCAL_ERROR="The interfaces on parent region number "// &
            & TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//" are not associated."        
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Parent region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("INTERFACE_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("INTERFACE_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE INTERFACE_USER_NUMBER_FIND
  
  !
  !================================================================================================================================
  !

  !>Finalises interfaces and deallocates all memory.
  SUBROUTINE INTERFACES_FINALISE(INTERFACES,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACES_TYPE), POINTER :: INTERFACES !<A pointer to the interfaces to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
     
    CALL ENTERS("INTERFACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACES)) THEN
      DO WHILE(INTERFACES%NUMBER_OF_INTERFACES>0)
        INTERFACE=>INTERFACES%INTERFACES(1)%PTR
        CALL INTERFACE_DESTROY(INTERFACE,ERR,ERROR,*999)
      ENDDO
      IF(ASSOCIATED(INTERFACES%INTERFACES)) DEALLOCATE(INTERFACES%INTERFACES)
      DEALLOCATE(INTERFACES)
    ENDIF
    
    CALL EXITS("INTERFACES_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACES_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACES_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises interfaces for a region.
  SUBROUTINE INTERFACES_INITIALISE(REGION,ERR,ERROR,*) 

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the interfaces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    CALL ENTERS("INTERFACES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%INTERFACES)) THEN
        LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has interfaces associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        ALLOCATE(REGION%INTERFACES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate region interfaces.",ERR,ERROR,*999)
        REGION%INTERFACES%PARENT_REGION=>REGION
        REGION%INTERFACES%NUMBER_OF_INTERFACES=0
        NULLIFY(REGION%INTERFACES%INTERFACES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACES_INITIALISE")
    RETURN
999 CALL INTERFACES_FINALISE(REGION%INTERFACES,ERR,ERROR,*998)
998 CALL ERRORS("INTERFACES_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACES_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface element connectivity.
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_INITIALISE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the mesh connectivity for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: InterfaceElementIdx,CoupledMeshIdx
     
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%MESH_CONNECTIVITY)) THEN
        IF(ALLOCATED(INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
          CALL FLAG_ERROR("Interface mesh element connectivity is already allocated.",ERR,ERROR,*999)
        ELSE
          IF(INTERFACE%NUMBER_OF_COUPLED_MESHES>0) THEN
            IF(INTERFACE%MESH_CONNECTIVITY%INTERFACE_MESH%NUMBER_OF_ELEMENTS>0) THEN
              INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS=INTERFACE%MESHES%MESHES(1)%PTR%NUMBER_OF_ELEMENTS
              INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES=INTERFACE%NUMBER_OF_COUPLED_MESHES
              IF (INTERFACE%SELF_CONTACT) THEN
                INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES=INTERFACE%NUMBER_OF_COUPLED_MESHES+1
              ENDIF
              ALLOCATE(INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS, &
                & INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface element connectivity.",ERR,ERROR,*999)
              DO InterfaceElementIdx=1,INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS
                DO CoupledMeshIdx=1,INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
                  INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)%COUPLED_MESH_ELEMENT_NUMBER=0
                ENDDO !CoupledMeshIdx
              ENDDO !InterfaceElementIdx
            ELSE
              CALL FLAG_ERROR("Interface coupled meshes are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface coupled meshes are not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface mesh connectivity is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_INITIALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface mesh connectivity.
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_POINTS_INITIALISE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the mesh connectivity for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: I, J, K
     
    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_POINTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%POINTS_CONNECTIVITY)) THEN
        IF(ALLOCATED(INTERFACE%POINTS_CONNECTIVITY%POINTS_CONNECTIVITY)) THEN
          CALL FLAG_ERROR("Interface points element connectivity is already allocated.",ERR,ERROR,*999)
        ELSE
          IF(INTERFACE%NUMBER_OF_COUPLED_MESHES<=0) CALL FLAG_ERROR("Interface coupled meshes are not associated.",ERR,ERROR,*999)
          !Hard coded the first mesh to be interface mesh
          IF(INTERFACE%MESHES%MESHES(1)%PTR%NUMBER_OF_ELEMENTS<=0) CALL FLAG_ERROR("Interface coupled meshes are not & 
            & associated.",ERR,ERROR,*999)
          INTERFACE%POINTS_CONNECTIVITY%NUMBER_OF_DATA_POINTS=INTERFACE%DATA_POINTS%NUMBER_OF_DATA_POINTS
          INTERFACE%POINTS_CONNECTIVITY%NUMBER_INT_DOM=INTERFACE%NUMBER_OF_COUPLED_MESHES
          ALLOCATE(INTERFACE%POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(INTERFACE%POINTS_CONNECTIVITY%NUMBER_OF_DATA_POINTS, &
            & INTERFACE%POINTS_CONNECTIVITY%NUMBER_INT_DOM))
          DO I = 1, INTERFACE%POINTS_CONNECTIVITY%NUMBER_OF_DATA_POINTS
            DO J = 1, INTERFACE%POINTS_CONNECTIVITY%NUMBER_INT_DOM
              INTERFACE%POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(I,J)%COUPLED_MESH_ELEMENT_NUMBER=0
            ENDDO! J
          ENDDO! I
        END IF
      ELSE
        CALL FLAG_ERROR("Interface points connectivity is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_POINTS_INITIALISE")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_POINTS_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_POINTS_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_POINTS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises an interface element connectivity and deallocates all memory.
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_FINALISE(INTERFACE_MESH_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE) :: INTERFACE_MESH_CONNECTIVITY !<The interface element connectivity to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: InterfaceElementIdx,CoupledMeshIdx
     
    CALL ENTERS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_FINALISE",ERR,ERROR,*999)

    DO InterfaceElementIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS  
      DO CoupledMeshIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
        IF(ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
          INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)%COUPLED_MESH_ELEMENT_NUMBER=0
          IF(ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)%XI)) THEN
            DEALLOCATE(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)%XI)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mesh connectivity element connectivity is being deallocated before allocation.", &
            & ERR,ERROR,*999)
        ENDIF
      ENDDO !InterfaceElementIdx
    ENDDO !CoupledMeshIdx

    DEALLOCATE(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)
    
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises an interface element connectivity and deallocates all memory.
  SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_POINTS_FINALISE(INTERFACE_POINTS_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE) :: INTERFACE_POINTS_CONNECTIVITY !<The interface element connectivity to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: I, J 
     
    CALL ENTERS("INTERFACE_POINTS_CONNECTIVITY_POINTS_FINALISE",ERR,ERROR,*999)
    !Deallocate data points connectivity 
    DO I = 1,INTERFACE_POINTS_CONNECTIVITY%NUMBER_OF_DATA_POINTS  !<Loop over all the interface data points
      DO J = 1,INTERFACE_POINTS_CONNECTIVITY%NUMBER_INT_DOM  !<Loop over all coupled domains 
        IF (ALLOCATED(INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY)) THEN
          INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(I,J)%COUPLED_MESH_ELEMENT_NUMBER=0
          INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(I,J)%COUPLED_MESH_CONTACT_NUMBER=0
          INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(I,J)%COUPLED_MESH_CONTACT_XI_NORMAL=0
          IF(ALLOCATED(INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(I,J)%XI)) THEN
            DEALLOCATE(INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(I,J)%XI)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mesh connectivity points connectivity is already deallocated.",ERR,ERROR,*999)
        ENDIF
      ENDDO ! END J - Loop over all coupled domains
    ENDDO ! END I - Loop over all the interface data points
    !Deallocate data points-coupled mesh element connectivity 
    DO I = 1,INTERFACE_POINTS_CONNECTIVITY%NUMBER_OF_ELEMENTS  !<Loop over all the interface elements
      DO J = 1,INTERFACE_POINTS_CONNECTIVITY%NUMBER_INT_DOM  !<Loop over all coupled domains 
        IF (ALLOCATED(INTERFACE_POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS)) THEN
          INTERFACE_POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS(I,J)%NUMBER_OF_COUPLED_MESH_ELEMENTS=0
          IF(ALLOCATED(INTERFACE_POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS(I,J)%ELEMENT_NUMBERS)) THEN
            DEALLOCATE(INTERFACE_POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS(I,J)%ELEMENT_NUMBERS)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mesh connectivity points connectivity is already deallocated.",ERR,ERROR,*999)
        ENDIF
      ENDDO ! END J - Loop over all coupled domains
    ENDDO ! END I - Loop over all the interface data points
    DEALLOCATE(INTERFACE_POINTS_CONNECTIVITY%POINTS_CONNECTIVITY)
    
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_POINTS_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_POINTS_CONNECTIVITY_POINTS_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_POINTS_CONNECTIVITY_POINTS_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_POINTS_CONNECTIVITY_POINTS_FINALISE

  !
  !================================================================================================================================
  !

END MODULE INTERFACE_ROUTINES

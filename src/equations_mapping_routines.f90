!> \file
!> \author Chris Bradley
!> \brief This module handles all equations mapping routines
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
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>This module handles all equations mapping routines.
MODULE EQUATIONS_MAPPING_ROUTINES

  USE BASE_ROUTINES
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES
 
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET
    MODULE PROCEDURE EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_1
    MODULE PROCEDURE EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_2
  END INTERFACE !EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET
  
  INTERFACE EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET
    MODULE PROCEDURE EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1
    MODULE PROCEDURE EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2
  END INTERFACE !EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET
  
  PUBLIC EQUATIONS_MAPPING_CREATE_FINISH,EQUATIONS_MAPPING_CREATE_START,EQUATIONS_MAPPING_DESTROY, &
    & EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET,EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET, &
    & EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET,EquationsMapping_LhsVariableTypeSet, &
    & EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFS_SET, &
    & EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET,EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET, &
    & EQUATIONS_MAPPING_RESIDUAL_COEFFS_SET,EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET, &
    & EQUATIONS_MAPPING_RESIDUAL_VARIABLES_NUMBER_SET, &
    & EQUATIONS_MAPPING_RHS_COEFF_SET,EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET, &
    & EQUATIONS_MAPPING_SOURCE_COEFF_SET,EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the equations/dofs mapping.
  SUBROUTINE EquationsMapping_Calculate(equationsMapping,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping !<A pointer to the equations mapping to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,dofIdx,matrixIdx,numberOfRows,numberOfGlobalRows,rowIdx, &
      & totalNumberOfRows,residualIdx,variableType
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: createValuesCache
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping
    TYPE(EquationsMapping_LhsType), POINTER :: lhsMapping
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: rhsMapping
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: sourceMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(FIELD_TYPE), POINTER :: dependentField,sourceField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable,lhsVariable,rhsVariable,sourceVariable
    TYPE(VARYING_STRING) :: localError

    CALL Enters("EquationsMapping_Calculate",err,error,*999)

    IF(ASSOCIATED(equationsMapping)) THEN
      equations=>equationsMapping%equations
      IF(ASSOCIATED(equations)) THEN
        createValuesCache=>equationsMapping%CREATE_VALUES_CACHE
        IF(ASSOCIATED(createValuesCache)) THEN
          equationsSet=>Equations%EQUATIONS_SET
          IF(ASSOCIATED(equationsSet)) THEN
            dependentField=>equationsSet%dependent%DEPENDENT_FIELD
            IF(ASSOCIATED(dependentField)) THEN              
              IF(createValuesCache%SOURCE_VARIABLE_TYPE/=0) THEN
                IF(ASSOCIATED(equationsSet%source)) THEN
                  sourceField=>equationsSet%source%SOURCE_FIELD
                  IF(.NOT.ASSOCIATED(sourceField)) CALL FlagError("Source field is not associated.",err,error,*999)
                ELSE
                  CALL FlagError("Equations set source is not associated.",err,error,*999)
                ENDIF
              ENDIF
              !Calculate the number of rows in the equations set
              NULLIFY(lhsVariable)
              CALL FIELD_VARIABLE_GET(dependentField,createValuesCache%lhsVariableType,lhsVariable,err,error,*999)
              IF(ASSOCIATED(lhsVariable%DOMAIN_MAPPING)) THEN
!!TODO: In general the rows will be dependent on the solution method used for the equations set. For now,
!! just do FEM problems.
                numberOfRows=lhsVariable%NUMBER_OF_DOFS
                totalNumberOfRows=lhsVariable%TOTAL_NUMBER_OF_DOFS
                equationsMapping%ROW_DOFS_MAPPING=>lhsVariable%DOMAIN_MAPPING
                numberOfGlobalRows=equationsMapping%ROW_DOFS_MAPPING%NUMBER_OF_GLOBAL
              ELSE
                CALL FlagError("LHS variable domain mapping is not associated.",err,error,*999)
              ENDIF
              !Check that the number of rows are consistent with the RHS vector if it exists
              IF(createValuesCache%RHS_VARIABLE_TYPE/=0) THEN
                NULLIFY(rhsVariable)
                CALL FIELD_VARIABLE_GET(dependentField,createValuesCache%RHS_VARIABLE_TYPE,rhsVariable,err,error,*999)
                IF(rhsVariable%NUMBER_OF_DOFS/=numberOfRows) THEN
                  localError="Invalid equations set up. The number of rows in the LHS of the equations set ("// &
                    & TRIM(NumberToVstring(numberOfRows,"*",err,error))// &
                    & ") does not match the number of rows in the RHS ("// &
                    & TRIM(NumberToVstring(rhsVariable%NUMBER_OF_DOFS,"*",err,error))//")."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                IF(rhsVariable%TOTAL_NUMBER_OF_DOFS/=totalNumberOfRows) THEN
                  localError="Invalid equations set up. The total number of rows in the LHS of the equations set ("// &
                    & TRIM(NumberToVstring(totalNumberOfRows,"*",err,error))// &
                    & ") does not match the total number of rows in the RHS ("// &
                    & TRIM(NumberToVstring(rhsVariable%TOTAL_NUMBER_OF_DOFS,"*",err,error))//")."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDIF
              equationsMapping%NUMBER_OF_ROWS=numberOfRows
              equationsMapping%TOTAL_NUMBER_OF_ROWS=totalNumberOfRows
              equationsMapping%NUMBER_OF_GLOBAL_ROWS=numberOfGlobalRows
              !Calculate dynamic mappings
              IF(createValuesCache%DYNAMIC_VARIABLE_TYPE/=0) THEN
                CALL EQUATIONS_MAPPING_DYNAMIC_MAPPING_INITIALISE(equationsMapping,err,error,*999)
                dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
                IF(ASSOCIATED(dynamicMapping)) THEN
                  dynamicMapping%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=createValuesCache%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                  dynamicMapping%STIFFNESS_MATRIX_NUMBER=createValuesCache%DYNAMIC_STIFFNESS_MATRIX_NUMBER
                  dynamicMapping%DAMPING_MATRIX_NUMBER=createValuesCache%DYNAMIC_DAMPING_MATRIX_NUMBER
                  dynamicMapping%MASS_MATRIX_NUMBER=createValuesCache%DYNAMIC_MASS_MATRIX_NUMBER
                  !Initialise the variable type maps
                  ALLOCATE(dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate equations mapping variable to equations map.",err,error,*999)
                  DO variableType=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                    CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE( &
                      & dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType),err,error,*999)
                    dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%VARIABLE_INDEX=variableType
                    dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%VARIABLE_TYPE=variableType
                    dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%Variable=>dependentField% &
                      & VARIABLE_TYPE_MAP(variableType)%ptr
                  ENDDO !variable_type
                  dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(createValuesCache%DYNAMIC_VARIABLE_TYPE)% &
                    & NUMBER_OF_EQUATIONS_MATRICES=createValuesCache%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                  IF(createValuesCache%RHS_VARIABLE_TYPE/=0) dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                    & createValuesCache%RHS_VARIABLE_TYPE)%NUMBER_OF_EQUATIONS_MATRICES=-1
                  !Allocate and initialise the variable to equations matrices maps
                  DO variableType=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                    dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
                    IF(ASSOCIATED(dependentVariable)) THEN
                      IF(dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES==-1) THEN
                        !???
                      ELSE IF(dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                        ALLOCATE(dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%EQUATIONS_MATRIX_NUMBERS( &
                          & dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES),STAT=err)
                        IF(err/=0) &
                          & CALL FlagError("Could not allocate variable to equations matrices maps equations matrix numbers.", &
                          & err,error,*999)
                        ALLOCATE(dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%DOF_TO_COLUMNS_MAPS( &
                          & dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES),STAT=err)
                        IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps dof to columns map.", &
                          & err,error,*999)                
                        dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%EQUATIONS_MATRIX_NUMBERS=0
                        IF(createValuesCache%DYNAMIC_VARIABLE_TYPE==variableType) THEN
                          IF(createValuesCache%DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                            dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%EQUATIONS_MATRIX_NUMBERS( &
                            & createValuesCache%DYNAMIC_STIFFNESS_MATRIX_NUMBER)=createValuesCache% &
                            & DYNAMIC_STIFFNESS_MATRIX_NUMBER
                          ENDIF
                          IF(createValuesCache%DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                            dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%EQUATIONS_MATRIX_NUMBERS( &
                            & createValuesCache%DYNAMIC_DAMPING_MATRIX_NUMBER)=createValuesCache% &
                            & DYNAMIC_DAMPING_MATRIX_NUMBER
                          ENDIF
                          IF(createValuesCache%DYNAMIC_MASS_MATRIX_NUMBER/=0) THEN
                            dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%EQUATIONS_MATRIX_NUMBERS( &
                            & createValuesCache%DYNAMIC_MASS_MATRIX_NUMBER)=createValuesCache% &
                            & DYNAMIC_MASS_MATRIX_NUMBER
                          ENDIF
                          DO matrixIdx=1,dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES
                            ALLOCATE(dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%DOF_TO_COLUMNS_MAPS( &
                              & matrixIdx)%COLUMN_DOF(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate variable dof to columns map column dof.", &
                              & err,error,*999)
                            DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                              !1-1 mapping for now
                              columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
                              dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%DOF_TO_COLUMNS_MAPS(matrixIdx)% &
                                & COLUMN_DOF(dofIdx)=columnIdx
                            ENDDO !dofIdx      
                          ENDDO !matrixIdx
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO !variableType
                  !Allocate and initialise the equations matrix to variable maps types
                  ALLOCATE(dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(dynamicMapping%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES), &
                    & STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate equations mapping equations matrix to variable maps.", &
                    & err,error,*999)
                  !Create the individual matrix maps and column maps
                  variableType=createValuesCache%DYNAMIC_VARIABLE_TYPE
                  dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
                  dynamicMapping%DYNAMIC_VARIABLE_TYPE=variableType
                  dynamicMapping%DYNAMIC_VARIABLE=>dependentVariable
                  DO matrixIdx=1,dynamicMapping%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                    CALL EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE(dynamicMapping% &
                      & EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx),err,error,*999)
                    dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%MATRIX_NUMBER=matrixIdx
                    dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%VARIABLE_TYPE=variableType
                    dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%variable=>dependentVariable
                    dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%NUMBER_OF_COLUMNS=dependentVariable% &
                      & DOMAIN_MAPPING%NUMBER_OF_GLOBAL
                    dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%MATRIX_COEFFICIENT=createValuesCache% &
                      & DYNAMIC_MATRIX_COEFFICIENTS(matrixIdx)
                    ALLOCATE(dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_TO_DOF_MAP( &
                      & dependentVariable%DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate equation matrix to variable map column to dof map.",&
                      & err,error,*999)
                    dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_TO_DOF_MAP=0
                    DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                      !1-1 mapping for now
                      columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
                      dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_TO_DOF_MAP(columnIdx)=dofIdx
                    ENDDO !dofIdx
                    dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_DOFS_MAPPING=>dependentVariable%DOMAIN_MAPPING
                  ENDDO !matrixIdx
                ELSE
                  CALL FlagError("Dynamic mapping is not associated.",err,error,*999)
                ENDIF
              ENDIF
              !Calculate linear mappings
              IF(createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN                  
                CALL EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE(equationsMapping,err,error,*999)
                linearMapping=>equationsMapping%LINEAR_MAPPING
                IF(ASSOCIATED(linearMapping)) THEN
                  linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  !Allocate and initialise the variable type maps
                  ALLOCATE(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate equations mapping variable to equations map.",err,error,*999)
                  DO variableType=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                    CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE( &
                      & linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType),err,error,*999)
                    linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%VARIABLE_INDEX=variableType
                    linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%VARIABLE_TYPE=variableType
                    linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%VARIABLE=>dependentField% &
                      & VARIABLE_TYPE_MAP(variableType)%ptr
                  ENDDO !
                  !Calculate the number of variable type maps and initialise
                  DO matrixIdx=1,linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    variableType=createValuesCache%LINEAR_MATRIX_VARIABLE_TYPES(matrixIdx)
                    linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES= &
                      & linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES+1
                  ENDDO !matrixIdx
                  IF(createValuesCache%RHS_VARIABLE_TYPE/=0) linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                    & createValuesCache%RHS_VARIABLE_TYPE)%NUMBER_OF_EQUATIONS_MATRICES=-1                    
                  linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES=0
                  !Allocate and initialise the variable to equations matrices maps
                  DO variableType=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                    dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
                    IF(ASSOCIATED(dependentVariable)) THEN
                      IF(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES==-1) THEN
                        linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES=linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES+1
                      ELSE IF(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                        ALLOCATE(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%EQUATIONS_MATRIX_NUMBERS( &
                          & linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES),STAT=err)
                        IF(err/=0) &
                          & CALL FlagError("Could not allocate variable to equations matrices maps equations matrix numbers.", &
                          & err,error,*999)
                        ALLOCATE(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType )%DOF_TO_COLUMNS_MAPS( &
                          & linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES),STAT=err)
                        IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps dof to columns map.", &
                          & err,error,*999)                
                        linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%EQUATIONS_MATRIX_NUMBERS=0
                        linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES=0
                        DO matrixIdx=1,linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          IF(createValuesCache%LINEAR_MATRIX_VARIABLE_TYPES(matrixIdx)==variableType) THEN
                            linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES= &
                              & linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES+1
                            linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%EQUATIONS_MATRIX_NUMBERS( &
                              & linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES)= &
                              & matrixIdx
                            ALLOCATE(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%DOF_TO_COLUMNS_MAPS( &
                              & linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                              & NUMBER_OF_EQUATIONS_MATRICES)%COLUMN_DOF(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate variable dof to columns map column dof.", &
                              & err,error,*999)
                            DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                              !1-1 mapping for now
                              columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
                              linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%DOF_TO_COLUMNS_MAPS( &
                                & linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                & NUMBER_OF_EQUATIONS_MATRICES)%COLUMN_DOF(dofIdx)=columnIdx
                            ENDDO !dofIdx
                          ENDIF
                        ENDDO !matrixIdx
                        linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES=linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES+1
                      ENDIF
                    ENDIF
                  ENDDO !variableType
                  !Allocate and initialise the variable types
                  ALLOCATE(linearMapping%LINEAR_MATRIX_VARIABLE_TYPES(linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate equations mapping matrix variable types.",err,error,*999)
                  linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES=0
                  DO variableType=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                    IF(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                      linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES=linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES+1
                      linearMapping%LINEAR_MATRIX_VARIABLE_TYPES(linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES)=variableType
                    ENDIF
                  ENDDO !variableType
                  !Allocate and initialise the equations matrix to variable maps types
                  ALLOCATE(linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES), &
                    & STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate equations mapping equations matrix to variable maps.", &
                    & err,error,*999)
                  !Create the individual matrix maps and column maps
                  DO matrixIdx=1,linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    variableType=createValuesCache%LINEAR_MATRIX_VARIABLE_TYPES(matrixIdx)
                    dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
                    CALL EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE(linearMapping% &
                      & EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx),err,error,*999)
                    linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%MATRIX_NUMBER=matrixIdx
                    linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%VARIABLE_TYPE=variableType
                    linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%variable=>dependentVariable
                    linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%NUMBER_OF_COLUMNS=dependentVariable% &
                      & DOMAIN_MAPPING%NUMBER_OF_GLOBAL
                    linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%MATRIX_COEFFICIENT=createValuesCache% &
                      & LINEAR_MATRIX_COEFFICIENTS(matrixIdx)
                    ALLOCATE(linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_TO_DOF_MAP( &
                      & dependentVariable%DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=err)                  
                    IF(err/=0) CALL FlagError("Could not allocate equation matrix to variable map column to dof map.",&
                      & err,error,*999)
                    linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_TO_DOF_MAP=0
                    DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                      !1-1 mapping for now
                      columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
                      linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_TO_DOF_MAP(columnIdx)=dofIdx
                    ENDDO !dofIdx
                    linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_DOFS_MAPPING=>dependentVariable%DOMAIN_MAPPING
                  ENDDO !matrixIdx
                ELSE
                  CALL FlagError("Linear mapping is not associated.",err,error,*999)
                ENDIF
              ENDIF
              !Calculate non-linear mappings
              IF(createValuesCache%NUMBER_OF_RESIDUAL_VARIABLES/=0) THEN
                CALL EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE(equationsMapping,err,error,*999)
                nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
                IF(ASSOCIATED(nonlinearMapping)) THEN
                  nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES=createValuesCache%NUMBER_OF_RESIDUAL_VARIABLES
                  ALLOCATE(nonlinearMapping%VAR_TO_JACOBIAN_MAP(nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate variable to Jacobian maps.",err,error,*999)
                  ALLOCATE(nonlinearMapping%JACOBIAN_TO_VAR_MAP(nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate Jacobian to variable maps.",err,error,*999)
                  ALLOCATE(nonlinearMapping%RESIDUAL_VARIABLES(nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate nonlinear mapping residual variables.",err,error,*999)
                  ALLOCATE(nonlinearMapping%RESIDUAL_COEFFICIENTS(nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate nonlinear mapping residual coefficients.",err,error,*999)
                  DO matrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                    CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE(nonlinearMapping% &
                      & VAR_TO_JACOBIAN_MAP(matrixIdx),err,error,*999)
                    nonlinearMapping%VAR_TO_JACOBIAN_MAP(matrixIdx)%JACOBIAN_NUMBER=matrixIdx
                    nonlinearMapping%VAR_TO_JACOBIAN_MAP(matrixIdx)%VARIABLE_TYPE= &
                      & createValuesCache%RESIDUAL_VARIABLE_TYPES(matrixIdx)
                    dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache% &
                      & RESIDUAL_VARIABLE_TYPES(matrixIdx))%ptr
                    nonlinearMapping%VAR_TO_JACOBIAN_MAP(matrixIdx)%variable=>dependentVariable
                    nonlinearMapping%RESIDUAL_VARIABLES(matrixIdx)%ptr=>dependentVariable
                    nonlinearMapping%RESIDUAL_COEFFICIENTS(matrixIdx)=createValuesCache%RESIDUAL_COEFFICIENTS(matrixIdx)
                    !Allocate and set dof to Jacobian columns map
                    ALLOCATE(nonlinearMapping%VAR_TO_JACOBIAN_MAP(matrixIdx)%DOF_TO_COLUMNS_MAP(dependentVariable% &
                      & TOTAL_NUMBER_OF_DOFS),STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate variable to Jacobian map dof to columns map.",err,error,*999)
                    DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                      !1-1 mapping for now
                      columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
                      nonlinearMapping%VAR_TO_JACOBIAN_MAP(matrixIdx)%DOF_TO_COLUMNS_MAP(dofIdx)=columnIdx
                    ENDDO !dofIdx
                    CALL EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE(nonlinearMapping% &
                      & JACOBIAN_TO_VAR_MAP(matrixIdx),err,error,*999)
                    nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%JACOBIAN_NUMBER=matrixIdx
                    nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%VARIABLE_TYPE=createValuesCache% &
                      & RESIDUAL_VARIABLE_TYPES(matrixIdx)
                    nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%variable=>dependentVariable
                    nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%NUMBER_OF_COLUMNS= &
                      & dependentVariable%DOMAIN_MAPPING%NUMBER_OF_GLOBAL
                    nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%JACOBIAN_COEFFICIENT= &
                      & createValuesCache%RESIDUAL_COEFFICIENTS(matrixIdx)
                    ALLOCATE(nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP( &
                      & dependentVariable%DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=err)
                    nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP=0
                    DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                      !1-1 mapping for now
                      columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
                      nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP(columnIdx)=dofIdx
                    ENDDO !dofIdx
                    nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%COLUMN_DOFS_MAPPING=>dependentVariable%DOMAIN_MAPPING
                  ENDDO !matrixIdx
                ELSE
                  CALL FlagError("Nonlinear mapping is not associated.",err,error,*999)
                ENDIF
              ENDIF
              !Calculate LHS mappings
              IF(createValuesCache%lhsVariableType==0) THEN
                CALL FlagError("Invalid equations setup. No LHS variable has been set.",err,error,*999)
              ELSE
                CALL EquationsMapping_LhsMappingInitialise(equationsMapping,err,error,*999)
                lhsMapping=>equationsMapping%lhsMapping
                IF(ASSOCIATED(lhsMapping)) THEN
                  lhsMapping%lhsVariableType=createValuesCache%lhsVariableType
                  dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%lhsVariableType)%ptr
                  lhsMapping%lhsVariable=>dependentVariable
                  lhsMapping%lhsVariableMapping=>dependentVariable%DOMAIN_MAPPING
                  !Allocate and set up the row mappings
                  ALLOCATE(lhsMapping%lhsDofToEquationsRowMap(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate LHS DOF to equations row map.",err,error,*999)
                  ALLOCATE(lhsMapping%equationsRowToLhsDofMap(totalNumberOfRows),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate equations row to LHS DOF map.",err,error,*999)
                  DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                    !1-1 mapping for now
                    rowIdx=dofIdx
                    lhsMapping%lhsDofToEquationsRowMap(dofIdx)=rowIdx
                  ENDDO !dofIdx
                  DO rowIdx=1,totalNumberOfRows
                    !1-1 mapping for now
                    dofIdx=rowIdx
                    lhsMapping%equationsRowToLhsDofMap(rowIdx)=dofIdx
                  ENDDO !rowIdx
                ELSE
                  CALL FlagError("LHS mapping is not associated.",err,error,*999)
                ENDIF
              ENDIF
              !Calculate RHS mappings
              IF(createValuesCache%RHS_VARIABLE_TYPE/=0) THEN                  
                CALL EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE(equationsMapping,err,error,*999)
                rhsMapping=>equationsMapping%RHS_MAPPING
                IF(ASSOCIATED(rhsMapping)) THEN
                  rhsMapping%RHS_VARIABLE_TYPE=createValuesCache%RHS_VARIABLE_TYPE
                  dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%RHS_VARIABLE_TYPE)%ptr
                  rhsMapping%RHS_VARIABLE=>dependentVariable
                  rhsMapping%RHS_VARIABLE_MAPPING=>dependentVariable%DOMAIN_MAPPING
                  rhsMapping%RHS_COEFFICIENT=createValuesCache%RHS_COEFFICIENT
                  !Allocate and set up the row mappings
                  ALLOCATE(rhsMapping%RHS_DOF_TO_EQUATIONS_ROW_MAP(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate RHS DOF to equations row map.",err,error,*999)
                  ALLOCATE(rhsMapping%EQUATIONS_ROW_TO_RHS_DOF_MAP(totalNumberOfRows),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate RHS equations row to DOF map.",err,error,*999)
                  DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                    !1-1 mapping for now
                    rowIdx=dofIdx
                    rhsMapping%RHS_DOF_TO_EQUATIONS_ROW_MAP(dofIdx)=rowIdx
                  ENDDO !dofIdx
                  DO rowIdx=1,totalNumberOfRows
                    !1-1 mapping for now
                    dofIdx=rowIdx
                    rhsMapping%EQUATIONS_ROW_TO_RHS_DOF_MAP(rowIdx)=dofIdx
                  ENDDO !rowIdx
                ELSE
                  CALL FlagError("RHS mapping is not associated.",err,error,*999)
                ENDIF
              ENDIF
              !Calcuate the source mappings
              IF(createValuesCache%SOURCE_VARIABLE_TYPE/=0) THEN                  
                CALL EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE(equationsMapping,err,error,*999)
                sourceMapping=>equationsMapping%SOURCE_MAPPING
                IF(ASSOCIATED(sourceMapping)) THEN
                  sourceMapping%SOURCE_VARIABLE_TYPE=createValuesCache%SOURCE_VARIABLE_TYPE
                  sourceVariable=>sourceField%VARIABLE_TYPE_MAP(createValuesCache%SOURCE_VARIABLE_TYPE)%ptr
                  sourceMapping%SOURCE_VARIABLE=>sourceVariable
                  sourceMapping%SOURCE_COEFFICIENT=createValuesCache%SOURCE_COEFFICIENT
                ELSE
                  CALL FlagError("Source mapping is not associated.",err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations equations set is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated.",err,error,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Equations mappings:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",equationsMapping%NUMBER_OF_ROWS,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of rows = ",equationsMapping%TOTAL_NUMBER_OF_ROWS, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global rows = ",equationsMapping%NUMBER_OF_GLOBAL_ROWS, &
        & err,error,*999)
      dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
      IF(ASSOCIATED(dynamicMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Dynamic mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dynamic equations matrices = ",dynamicMapping% &
          & NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic stiffness matrix number = ",dynamicMapping% &
          & STIFFNESS_MATRIX_NUMBER,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic damping matrix number = ",dynamicMapping% &
          & DAMPING_MATRIX_NUMBER,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic mass matrix number = ",dynamicMapping% &
          & MASS_MATRIX_NUMBER,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic variable type = ",dynamicMapping% &
          & DYNAMIC_VARIABLE_TYPE,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",err,error,*999)
        DO variableType=1,FIELD_NUMBER_OF_VARIABLE_TYPES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type : ",variableType,err,error,*999)
          IF(ASSOCIATED(dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%variable)) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",dynamicMapping% &
              & VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%variable%TOTAL_NUMBER_OF_DOFS,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations matrices = ",dynamicMapping% &
              & VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES,err,error,*999)
            IF(dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                & NUMBER_OF_EQUATIONS_MATRICES,4,4,dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                & EQUATIONS_MATRIX_NUMBERS,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))',err,error,*999) 
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      DOF to column maps :",err,error,*999)
              DO matrixIdx=1,dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix number : ",matrixIdx,err,error,*999)
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                  & variableType)%VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5,dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                  & variableType)%DOF_TO_COLUMNS_MAPS(matrixIdx)%COLUMN_DOF, &
                  & '("        Column numbers :",5(X,I13))','(24X,5(X,I13))',err,error,*999) 
              ENDDO !matrixIdx
            ENDIF
          ENDIF
        ENDDO !variableType
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",err,error,*999)
        DO matrixIdx=1,dynamicMapping%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",dynamicMapping% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%VARIABLE_TYPE,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of columns = ",dynamicMapping% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%NUMBER_OF_COLUMNS,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",dynamicMapping% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%MATRIX_COEFFICIENT,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)% &
            & NUMBER_OF_COLUMNS,5,5,dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_TO_DOF_MAP, &
            & '("        Column to DOF maps :",5(X,I13))','(28X,5(X,I13))',err,error,*999) 
        ENDDO !matrixIdx
      ENDIF
      linearMapping=>equationsMapping%LINEAR_MAPPING
      IF(ASSOCIATED(linearMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Linear mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of linear equations matrices = ",linearMapping% &
          & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of linear matrix variables = ",linearMapping% &
          & NUMBER_OF_LINEAR_MATRIX_VARIABLES,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%NUMBER_OF_LINEAR_MATRIX_VARIABLES,4,4, &
          & linearMapping%LINEAR_MATRIX_VARIABLE_TYPES,'("    Linear matrix variable types :",4(X,I12))','(27X,4(X,I12))', &
          & err,error,*999) 
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",err,error,*999)
        DO variableType=1,FIELD_NUMBER_OF_VARIABLE_TYPES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type : ",variableType,err,error,*999)
          IF(ASSOCIATED(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%variable)) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",linearMapping% &
              & VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%VARIABLE%TOTAL_NUMBER_OF_DOFS,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations matrices = ",linearMapping% &
              & VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES,err,error,*999)
            IF(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                & NUMBER_OF_EQUATIONS_MATRICES,4,4,linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                & EQUATIONS_MATRIX_NUMBERS,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))',err,error,*999) 
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      DOF to column maps :",err,error,*999)
              DO matrixIdx=1,linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix number : ",matrixIdx,err,error,*999)
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                  & variableType)%VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5,linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                  & variableType)%DOF_TO_COLUMNS_MAPS(matrixIdx)%COLUMN_DOF, &
                  & '("        Column numbers :",5(X,I13))','(24X,5(X,I13))',err,error,*999) 
              ENDDO !matrixIdx
            ENDIF
          ENDIF
        ENDDO !variableType
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",err,error,*999)
        DO matrixIdx=1,linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",linearMapping% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%VARIABLE_TYPE,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of columns = ",linearMapping% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%NUMBER_OF_COLUMNS,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",linearMapping% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%MATRIX_COEFFICIENT,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)% &
            & NUMBER_OF_COLUMNS,5,5,linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%COLUMN_TO_DOF_MAP, &
            & '("        Column to DOF maps :",5(X,I13))','(28X,5(X,I13))',err,error,*999) 
        ENDDO !matrixIdx
      ENDIF
      nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
      IF(ASSOCIATED(nonlinearMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Nonlinear mappings:",err,error,*999)
        DO residualIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Residual number : ",residualIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Residual variable type = ",nonlinearMapping% &
            & JACOBIAN_TO_VAR_MAP(residualIdx)%VARIABLE_TYPE,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of residual DOFs = ",nonlinearMapping% &
            & JACOBIAN_TO_VAR_MAP(residualIdx)%VARIABLE%TOTAL_NUMBER_OF_DOFS,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Residual coefficient = ",nonlinearMapping%RESIDUAL_COEFFICIENTS( &
            & residualIdx),err,error,*999)
       ENDDO
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Jacobian mappings:",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to Jacobian mappings:",err,error,*999)
        DO matrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian variable type = ",nonlinearMapping% &
            & VAR_TO_JACOBIAN_MAP(matrixIdx)%VARIABLE_TYPE,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of Jacobain DOFs = ",nonlinearMapping% &
            & VAR_TO_JACOBIAN_MAP(matrixIdx)%VARIABLE%TOTAL_NUMBER_OF_DOFS,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nonlinearMapping%VAR_TO_JACOBIAN_MAP(matrixIdx)%variable% &
            & TOTAL_NUMBER_OF_DOFS,5,5,nonlinearMapping%VAR_TO_JACOBIAN_MAP(matrixIdx)%DOF_TO_COLUMNS_MAP, &
            & '("      DOF to column map :",5(X,I13))','(26X,5(X,I13))',err,error,*999) 
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Jacobian to variable mappings:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian variable type = ",nonlinearMapping% &
            & JACOBIAN_TO_VAR_MAP(matrixIdx)%VARIABLE_TYPE,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",nonlinearMapping% &
            & JACOBIAN_TO_VAR_MAP(matrixIdx)%NUMBER_OF_COLUMNS,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian coefficient = ",nonlinearMapping% &
            & JACOBIAN_TO_VAR_MAP(matrixIdx)%JACOBIAN_COEFFICIENT,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%NUMBER_OF_COLUMNS, &
            & 5,5,nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP, &
            & '("      Column to DOF map :",5(X,I13))','(26X,5(X,I13))',err,error,*999) 
        ENDDO
      ENDIF
      lhsMapping=>equationsMapping%lhsMapping
      IF(ASSOCIATED(lhsMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  LHS mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    LHS variable type = ",lhsMapping%lhsVariableType,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of LHS DOFs = ",lhsMapping%lhsVariable% &
          & TOTAL_NUMBER_OF_DOFS,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,lhsMapping%lhsVariable%TOTAL_NUMBER_OF_DOFS,5,5, &
          & lhsMapping%lhsDofToEquationsRowMap,'("    DOF to row mappings :",5(X,I13))','(25X,5(X,I13))',err,error,*999) 
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsMapping%TOTAL_NUMBER_OF_ROWS,5,5, &
          & lhsMapping%equationsRowToLhsDofMap,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))',err,error,*999) 
      ENDIF
      rhsMapping=>equationsMapping%RHS_MAPPING
      IF(ASSOCIATED(rhsMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  RHS mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    RHS variable type = ",rhsMapping%RHS_VARIABLE_TYPE,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of RHS DOFs = ",rhsMapping%RHS_VARIABLE% &
          & TOTAL_NUMBER_OF_DOFS,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    RHS coefficient = ",rhsMapping%RHS_COEFFICIENT,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,rhsMapping%RHS_VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5, &
          & rhsMapping%RHS_DOF_TO_EQUATIONS_ROW_MAP,'("    DOF to row mappings :",5(X,I13))','(25X,5(X,I13))',err,error,*999) 
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsMapping%TOTAL_NUMBER_OF_ROWS,5,5, &
          & rhsMapping%EQUATIONS_ROW_TO_RHS_DOF_MAP,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))',err,error,*999) 
      ENDIF
      sourceMapping=>equationsMapping%SOURCE_MAPPING
      IF(ASSOCIATED(sourceMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Source mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Source variable type = ",sourceMapping%SOURCE_VARIABLE_TYPE, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of source DOFs = ",sourceMapping%SOURCE_VARIABLE% &
          & TOTAL_NUMBER_OF_DOFS,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Source coefficient = ",sourceMapping%SOURCE_COEFFICIENT,err,error,*999)
      ENDIF
    ENDIF
       
    CALL Exits("EquationsMapping_Calculate")
    RETURN
999 CALL Errors("EquationsMapping_Calculate",err,error)
    CALL Exits("EquationsMapping_Calculate")
    RETURN 1
  END SUBROUTINE EquationsMapping_Calculate

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equations mapping
  SUBROUTINE EQUATIONS_MAPPING_CREATE_FINISH(equationsMapping,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping !<A pointer to the equations mapping
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx
    LOGICAL :: lhsIsUsed,isResidualType
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: createValuesCache
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(VARYING_STRING) :: localError

    CALL Enters("EQUATIONS_MAPPING_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(equationsMapping)) THEN
      IF(equationsMapping%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has already been finished.",err,error,*999)
      ELSE
        createValuesCache=>equationsMapping%CREATE_VALUES_CACHE
        IF(ASSOCIATED(createValuesCache)) THEN
          equations=>equationsMapping%equations
          IF(ASSOCIATED(equations)) THEN
            equationsSet=>equations%EQUATIONS_SET
            IF(ASSOCIATED(equationsSet)) THEN
              dependentField=>equationsSet%dependent%DEPENDENT_FIELD
              IF(ASSOCIATED(dependentField)) THEN              
                !Check that all the variables have been mapped properly
                SELECT CASE(equations%TIME_DEPENDENCE)
                CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
                  SELECT CASE(equations%linearity)
                  CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                    !Static linear equations set
                    IF(createValuesCache%RHS_VARIABLE_TYPE==0.AND.createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES==0) &
                      & CALL FlagError("Invalid equations mapping. The RHS variable type must be set if there are no "// &
                      & "linear matrices.",err,error,*999)
                    IF(createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>=1) THEN
                      !Use first listed nonlinear variable as LHS default
                      dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%LINEAR_MATRIX_VARIABLE_TYPES(1))%ptr
                    ELSE
                      CALL FlagError("The number of linear equations matrices must be at least one for a linear equations set.", &
                        & err,error,*999)
                    ENDIF
                  CASE(EQUATIONS_NONLINEAR)
                    !Static nonlinear equations set
                    DO matrixIdx=1,createValuesCache%NUMBER_OF_RESIDUAL_VARIABLES
                      IF(createValuesCache%RESIDUAL_VARIABLE_TYPES(matrixIdx)==0) THEN
                        localError="Invalid equations mapping. The residual variable type is not set for Jacobian number "// &
                          & TRIM(NumberToVstring(matrixIdx,"*",err,error))//"."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                    ENDDO
                    IF(createValuesCache%RHS_VARIABLE_TYPE==0.AND.createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES==0) &
                      & CALL FlagError("Invalid equations mapping. The RHS variable type must be set if there are no "// &
                      & "linear matrices.",err,error,*999)
                    !Use first listed nonlinear variable as LHS default
                    IF(createValuesCache%NUMBER_OF_RESIDUAL_VARIABLES>=1) THEN
                      dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%RESIDUAL_VARIABLE_TYPES(1))%ptr
                    ELSE
                      CALL FlagError("The number of Jacobian matrices must be at least one for a nonlinear equations set.", &
                        & err,error,*999)
                    ENDIF
                  CASE DEFAULT
                    localError="The equations linearity type of "//TRIM(NumberToVstring(equations%linearity,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
                  SELECT CASE(equations%linearity)
                  CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                    !Dynamic linear equations set
                    IF(createValuesCache%DYNAMIC_VARIABLE_TYPE==0) CALL FlagError("Invalid equations mapping. "// &
                      & "The dynamic variable type must be set for dynamic equations.", err,error,*999)
                    IF(createValuesCache%RHS_VARIABLE_TYPE==0.AND.createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES==0) &
                      & CALL FlagError("Invalid equations mapping. The RHS variable type must be set if there are no "// &
                      & "linear matrices.",err,error,*999)
                    !Use dynamic variable as LHS default
                    dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%DYNAMIC_VARIABLE_TYPE)%ptr
                  CASE(EQUATIONS_NONLINEAR)
                    !Dynamic nonlinear equations set
                    IF(createValuesCache%DYNAMIC_VARIABLE_TYPE==0) CALL FlagError("Invalid equations mapping. "// &
                      & "The dynamic variable type must be set for dynamic equations.", err,error,*999)
                    IF(createValuesCache%RHS_VARIABLE_TYPE==0.AND.createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES==0) &
                      & CALL FlagError("Invalid equations mapping. The RHS variable type must be set if there are no "// &
                      & "linear matrices.",err,error,*999)
                    isResidualType=.FALSE.
                    DO residualIdx=1,createValuesCache%NUMBER_OF_RESIDUAL_VARIABLES
                      IF(createValuesCache%RESIDUAL_VARIABLE_TYPES(residualIdx)==0) THEN
                        localError="Invalid equations mapping. The residual variable type is not set for residual number "// &
                          & TRIM(NumberToVstring(residualIdx,"*",err,error))//"."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                      isResidualType= &
                        & createValuesCache%RESIDUAL_VARIABLE_TYPES(residualIdx)==createValuesCache%DYNAMIC_VARIABLE_TYPE
                    ENDDO !residualIdx
                    IF(.NOT.isResidualType) CALL FlagError("Invalid equations mapping. "// &
                      & "The residual variable type must correspond to the dynamic variable "// &
                      & "type for nonlinear dynamic equations.", err,error,*999)
                    !Use first nonlinear/dynamic variable as LHS default
                    dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%RESIDUAL_VARIABLE_TYPES(1))%ptr
                  CASE DEFAULT
                    localError="The equations linearity type of "//TRIM(NumberToVstring(equations%linearity,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(EQUATIONS_TIME_STEPPING)
                  !Time stepping DAE equations set
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The equations time dependence type of "// &
                    & TRIM(NumberToVstring(equations%TIME_DEPENDENCE,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                !Check the linear matrices variable types
                DO matrixIdx=1,createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  IF(createValuesCache%LINEAR_MATRIX_VARIABLE_TYPES(matrixIdx)==0) THEN
                    localError="Invalid equations mapping. The linear matrix variable type is not set for linear matrix number "//&
                      & TRIM(NumberToVstring(matrixIdx,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDDO !matrix_idx
                !Check that the LHS has been set, if not default appropriately
                IF(createValuesCache%lhsVariableType/=0) THEN
                  !Check that LHS variable is used
                  lhsIsUsed=createValuesCache%lhsVariableType==createValuesCache%DYNAMIC_VARIABLE_TYPE
                  IF(.NOT.lhsIsUsed) THEN
                    DO matrixIdx=1,createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                      IF(createValuesCache%lhsVariableType==createValuesCache%LINEAR_MATRIX_VARIABLE_TYPES(matrixIdx)) THEN
                        lhsIsUsed=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !matrixIdx
                    IF(.NOT.lhsIsUsed) THEN
                      DO residualIdx=1,createValuesCache%NUMBER_OF_RESIDUAL_VARIABLES
                        IF(createValuesCache%lhsVariableType==createValuesCache%RESIDUAL_VARIABLE_TYPES(residualIdx)) THEN
                          lhsIsUsed=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO !residualIdx
                    ENDIF
                  ENDIF
                  IF(.NOT.lhsIsUsed) CALL FlagError("Invalid equations mapping. The dependent variable mapped to the LHS "// &
                    & "is not used on the LHS.",err,error,*999)
                ELSE
                  IF(ASSOCIATED(dependentVariable)) THEN
                    createValuesCache%lhsVariableType=dependentVariable%VARIABLE_TYPE
                  ELSE
                    CALL FlagError("The default dependent variable mapped to the LHS is not associated.",err,error,*999)
                  ENDIF
                ENDIF
                !Now calculate the equations mapping and clean up
                CALL EquationsMapping_Calculate(equationsMapping,err,error,*999)
                CALL EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE(equationsMapping%CREATE_VALUES_CACHE,err,error,*999)
                equationsMapping%EQUATIONS_MAPPING_FINISHED=.TRUE.
              ELSE
                CALL FlagError("The equations set dependent field is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations equations set is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated.",err,error,*999)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_CREATE_FINISH")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_CREATE_FINISH",err,error)
    CALL Exits("EQUATIONS_MAPPING_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equations mapping for a equations set equations
  SUBROUTINE EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equation to create the equations mapping from.
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<On return, a pointer to the equations mapping. This must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL Enters("EQUATIONS_MAPPING_CREATE_START",err,error,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          CALL FlagError("Equations mapping is already associated.",err,error,*999)
        ELSE
          NULLIFY(EQUATIONS_MAPPING)
          CALL EQUATIONS_MAPPING_INITIALISE(EQUATIONS,err,error,*999)
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        ENDIF
      ELSE
        CALL FlagError("Equations has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",err,error,*999)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_CREATE_START")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_CREATE_START",err,error)
    CALL Exits("EQUATIONS_MAPPING_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises an equations mapping create values cache and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE",err,error,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS)
      IF(ALLOCATED(CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS)
      IF(ALLOCATED(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENTS)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations mapping create values cache 
  SUBROUTINE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE(EQUATIONS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,matrix_idx2,VARIABLE_NUMBER
    LOGICAL :: IS_RESIDUAL_TYPE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FlagError("Equations mapping create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
          IF(ASSOCIATED(EQUATIONS_SET)) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              !Allocate and initialise the create values cache
              ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE,STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate equations mapping create values cache.",err,error,*999)
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%lhsVariableType=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT=1.0_DP
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE=0
              EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_COEFFICIENT=1.0_DP
              !Set the default equations mapping in the create values cache
              !First calculate how many linear and dynamic matrices we have and set the variable types for the dynamic, residual
              !and RHS variables
              IF(DEPENDENT_FIELD%NUMBER_OF_VARIABLES==1) THEN
                SELECT CASE(EQUATIONS%LINEARITY)
                CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                  CALL FlagError("Dependent field only has one variable which cannot be mapped to both an equations matrix "// &
                    & "and rhs vector.",err,error,*999)
                CASE(EQUATIONS_NONLINEAR)
                  CALL FlagError("Dependent field only has one variable which cannot be mapped to both the residual "// &
                    & "and rhs vector.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations linearity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",err,error))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE IF(DEPENDENT_FIELD%NUMBER_OF_VARIABLES>1) THEN
                SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
                CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
                  SELECT CASE(EQUATIONS%LINEARITY)
                  CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=DEPENDENT_FIELD%NUMBER_OF_VARIABLES-1
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%PTR)) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=DEPENDENT_FIELD% &
                        & VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                    ELSE
                      CALL FlagError("Not implemented.",err,error,*999)
                    ENDIF
                  CASE(EQUATIONS_NONLINEAR)
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES=1
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%PTR)) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=DEPENDENT_FIELD% &
                        & VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                    ELSE
                      CALL FlagError("Not implemented.",err,error,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The equations linearity type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",err,error))//" is invalid."
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  END SELECT
                CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
                  SELECT CASE(EQUATIONS%LINEARITY)
                  CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                    IF(EQUATIONS%TIME_DEPENDENCE==EQUATIONS_FIRST_ORDER_DYNAMIC) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=2
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=1
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=2
                    ELSE
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=3
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=1
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=2
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER=3
                    ENDIF
                    !EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=DEPENDENT_FIELD%NUMBER_OF_VARIABLES-2
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR)) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE=DEPENDENT_FIELD% &
                        & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                    ELSE
                      CALL FlagError("Not implemented.",err,error,*999)
                    ENDIF
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%PTR)) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=DEPENDENT_FIELD% &
                        & VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                    ELSE
                      CALL FlagError("Not implemented.",err,error,*999)
                    ENDIF
                  CASE(EQUATIONS_NONLINEAR)
! SEBK 19/08/2009 not sure about mapping here
!|
                    IF(EQUATIONS%TIME_DEPENDENCE==EQUATIONS_FIRST_ORDER_DYNAMIC) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=2
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=1
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=2
                    ELSE
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=3
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=1
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=2
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER=3
                    ENDIF
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
                    EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES=1
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR)) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE=DEPENDENT_FIELD% &
                        & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                    ELSE
                      CALL FlagError("Not implemented.",err,error,*999)
                    ENDIF
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%PTR)) THEN
                      EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=DEPENDENT_FIELD% &
                        & VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                    ELSE
                      CALL FlagError("Not implemented.",err,error,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The equations linearity type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",err,error))//" is invalid."
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  END SELECT
!|
! SEBK 19/08/2009 not sure about mapping here
                CASE DEFAULT
                  LOCAL_ERROR="The equations time dependence type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",err,error))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE
                LOCAL_ERROR="The number of dependent field variables of "// &
                  & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD%NUMBER_OF_VARIABLES,"*",err,error))//" is invalid."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              ENDIF
              !Allocate the dynamic matrix coefficients and set their values
              IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES>0) THEN
                ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(EQUATIONS_MAPPING% &
                  & CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) &
                  & CALL FlagError("Could not allocate equations mapping create values cache dynamic matrix coefficients.", &
                  & err,error,*999)
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS=1.0_DP !Equations matrices are added by default
              ENDIF
              !Allocate the residual variable types
              IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES>0) THEN
                ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES(EQUATIONS_MAPPING% &
                  & CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations mapping create values cache residual variable types.", &
                    & err,error,*999)
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES=0
                ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENTS(EQUATIONS_MAPPING% &
                  & CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations mapping create values cache residual coefficients.", &
                    & err,error,*999)
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENTS=1.0_DP
                DO matrix_idx=1,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES
                  VARIABLE_NUMBER=1
                  DO WHILE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES(matrix_idx)==0.AND. &
                    & VARIABLE_NUMBER<=FIELD_NUMBER_OF_VARIABLE_TYPES)
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_NUMBER)%PTR)) THEN
                      IF(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_NUMBER)%PTR%VARIABLE_TYPE/= &
                        & EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE) THEN
                        EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES(matrix_idx)= &
                          & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_NUMBER)%PTR%VARIABLE_TYPE
                      ENDIF
                    ENDIF
                    VARIABLE_NUMBER=VARIABLE_NUMBER+1
                  ENDDO
                ENDDO !matrix_idx
                IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES(EQUATIONS_MAPPING% &
                  & CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES)==0) &
                  & CALL FlagError("Invalid setup. All Jacobian matrices do not have a mapped dependent field variable.", &
                  & err,error,*999)
              ENDIF
              !Allocate the linear matrix variable types and linear matrix coefficients and set their values
              IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN
                ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING% &
                  & CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) CALL  & 
                  & FLAG_ERROR("Could not allocate equations mapping create values cache linear matrix variable types.", &
                  & err,error,*999)
                ALLOCATE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(EQUATIONS_MAPPING% &
                  & CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations mapping create values cache linear matrix coefficients.", &
                  & err,error,*999)
                !Set up the matrices variable types
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES=0
                VARIABLE_NUMBER=1
                DO matrix_idx=1,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  DO WHILE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)==0.AND. &
                    & VARIABLE_NUMBER<=FIELD_NUMBER_OF_VARIABLE_TYPES)
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_NUMBER)%PTR)) THEN
                      IF(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_NUMBER)%PTR%VARIABLE_TYPE/= &
                        & EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE) THEN
                        IS_RESIDUAL_TYPE=.FALSE.
                        DO matrix_idx2=1,EQUATIONS_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES
                          IF(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_NUMBER)%PTR%VARIABLE_TYPE== &
                            & EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES(matrix_idx2)) THEN
                            IS_RESIDUAL_TYPE=.TRUE.
                          ENDIF
                        ENDDO
                        IF(IS_RESIDUAL_TYPE.EQV..FALSE.) THEN
                          EQUATIONS_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)= &
                            & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_NUMBER)%PTR%VARIABLE_TYPE
                        ENDIF
                      ENDIF
                    ENDIF
                    VARIABLE_NUMBER=VARIABLE_NUMBER+1
                  ENDDO
                ENDDO !matrix_idx
                IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(EQUATIONS_MAPPING% &
                  & CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES)==0) &
                  & CALL FlagError("Invalid setup. All linear matrices do not have a mapped dependent field variable.", &
                  & err,error,*999)
                EQUATIONS_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS=1.0_DP !Equations matrices are added by default
              ENDIF
            ELSE
              CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("The equations equations set is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FlagError("The equations mapping equations is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL Errors("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Destroy an equations mapping.
  SUBROUTINE EQUATIONS_MAPPING_DESTROY(EQUATIONS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer the equations mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_DESTROY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      CALL EQUATIONS_MAPPING_FINALISE(EQUATIONS_MAPPING,err,error,*999)
    ELSE
      CALL FlagError("Equations mapping is not associated.",err,error,*999)
    ENDIF
        
    CALL Exits("EQUATIONS_MAPPING_DESTROY")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_DESTROY",err,error)    
    CALL Exits("EQUATIONS_MAPPING_DESTROY")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping dynamic  mapping and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MAPPING_FINALISE(DYNAMIC_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING !<A pointer to the dynamic mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,variable_type
 
    CALL Enters("EQUATIONS_MAPPING_DYNAMIC_MAPPING_FINALISE",err,error,*999)

    IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
      IF(ALLOCATED(DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS)) THEN
        DO variable_type=1,SIZE(DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS,1)
          CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE(DYNAMIC_MAPPING% &
            & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type),err,error,*999)
        ENDDO !variable_type
        DEALLOCATE(DYNAMIC_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS)        
      ENDIF
      IF(ALLOCATED(DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS)) THEN
        DO matrix_idx=1,SIZE(DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS,1)
          CALL EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE(DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx), &
            & err,error,*999)
        ENDDO !matrix_idx
        DEALLOCATE(DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS)
      ENDIF
      DEALLOCATE(DYNAMIC_MAPPING)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MAPPING_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_DYNAMIC_MAPPING_FINALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MAPPING_FINALISE")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping dynamic mapping
  SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MAPPING_INITIALISE(EQUATIONS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the dynamic mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL Enters("EQUATIONS_MAPPING_DYNAMIC_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%DYNAMIC_MAPPING)) THEN
        CALL FlagError("Equations mapping dynamic mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%DYNAMIC_MAPPING,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equations mapping dynamic mapping.",err,error,*999)
        EQUATIONS_MAPPING%DYNAMIC_MAPPING%EQUATIONS_MAPPING=>EQUATIONS_MAPPING
        EQUATIONS_MAPPING%DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=0
        EQUATIONS_MAPPING%DYNAMIC_MAPPING%STIFFNESS_MATRIX_NUMBER=0
        EQUATIONS_MAPPING%DYNAMIC_MAPPING%DAMPING_MATRIX_NUMBER=0
        EQUATIONS_MAPPING%DYNAMIC_MAPPING%MASS_MATRIX_NUMBER=0
        EQUATIONS_MAPPING%DYNAMIC_MAPPING%DYNAMIC_VARIABLE_TYPE=0
        NULLIFY(EQUATIONS_MAPPING%DYNAMIC_MAPPING%DYNAMIC_VARIABLE)
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_DYNAMIC_MAPPING_FINALISE(EQUATIONS_MAPPING%DYNAMIC_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL Errors("EQUATIONS_MAPPING_DYNAMIC_MAPPING_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MAPPING_INITIALISE")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrices involved in dynamic equations mapping
  SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_ALL(EQUATIONS_MAPPING,MASS_MATRIX,DAMPING_MATRIX,STIFFNESS_MATRIX, &
    & err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set the atrices for
    LOGICAL, INTENT(IN) :: MASS_MATRIX !<Is .TRUE. if the mass matrix is in the equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: DAMPING_MATRIX !<Is .TRUE. if the damping matrix is in the equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: STIFFNESS_MATRIX !<Is .TRUE. if the stiffness matrix is in the equations mapping, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_DYNAMIC_DAMPING_MATRIX_NUMBER,NEW_DYNAMIC_MASS_MATRIX_NUMBER,NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER, &
      & NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
    REAL(DP), ALLOCATABLE :: OLD_DYNAMIC_MATRIX_COEFFICIENTS(:)
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_ALL",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
        IF(ASSOCIATED(EQUATIONS_SET)) THEN
          CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
          IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
            SELECT CASE(EQUATIONS%LINEARITY)
            CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
              NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=0
              NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER=0
              NEW_DYNAMIC_DAMPING_MATRIX_NUMBER=0
              NEW_DYNAMIC_MASS_MATRIX_NUMBER=0
              IF(STIFFNESS_MATRIX) THEN
                NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES+1
                NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              ENDIF
              IF(DAMPING_MATRIX) THEN
                NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES+1
                NEW_DYNAMIC_DAMPING_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              ENDIF
              IF(MASS_MATRIX) THEN
                NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES+1
                NEW_DYNAMIC_MASS_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              ENDIF
              IF(NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES>0) THEN
                ALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate old dynamic matrix coefficients.",err,error,*999)
                OLD_DYNAMIC_MATRIX_COEFFICIENTS=CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS
                DEALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS)
                ALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate dynamic matrix coefficients.",err,error,*999)
                IF(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                IF(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                IF(NEW_DYNAMIC_MASS_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_MASS_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_MASS_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER
                CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=NEW_DYNAMIC_DAMPING_MATRIX_NUMBER
                CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER=NEW_DYNAMIC_MASS_MATRIX_NUMBER
                IF(ALLOCATED(OLD_DYNAMIC_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS)
              ELSE
                CALL FlagError("Invalid dynamic matrices set up. There are no dynamic equations matrices.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_NONLINEAR)
! SEBK 19/08/2009 not sure about mapping here
!|
              NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=0
              NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER=0
              NEW_DYNAMIC_DAMPING_MATRIX_NUMBER=0
              NEW_DYNAMIC_MASS_MATRIX_NUMBER=0
              IF(STIFFNESS_MATRIX) THEN
                NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES+1
                NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              ENDIF
              IF(DAMPING_MATRIX) THEN
                NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES+1
                NEW_DYNAMIC_DAMPING_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              ENDIF
              IF(MASS_MATRIX) THEN
                NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES+1
                NEW_DYNAMIC_MASS_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              ENDIF
              IF(NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES>0) THEN
                ALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate old dynamic matrix coefficients.",err,error,*999)
                OLD_DYNAMIC_MATRIX_COEFFICIENTS=CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS
                DEALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS)
                ALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate dynamic matrix coefficients.",err,error,*999)
                IF(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                IF(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                IF(NEW_DYNAMIC_MASS_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_MASS_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_MASS_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER
                CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=NEW_DYNAMIC_DAMPING_MATRIX_NUMBER
                CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER=NEW_DYNAMIC_MASS_MATRIX_NUMBER
                IF(ALLOCATED(OLD_DYNAMIC_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS)
              ELSE
                CALL FlagError("Invalid dynamic matrices set up. There are no dynamic equations matrices.",err,error,*999)
              ENDIF
!|
! SEBK 19/08/2009 not sure about mapping here
            CASE DEFAULT
              LOCAL_ERROR="The equations linearity type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Equations mapping create values cache is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations equations set is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_ALL_ORDER")
    RETURN
999 IF(ALLOCATED(OLD_DYNAMIC_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS)    
    CALL Errors("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_ALL",err,error)
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_ALL")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_ALL

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrices involved in a first order dynamic equations mapping
  SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_1(EQUATIONS_MAPPING,DAMPING_MATRIX,STIFFNESS_MATRIX,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set the first order matrices for
    LOGICAL, INTENT(IN) :: DAMPING_MATRIX !<Is .TRUE. if the damping matrix is in the equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: STIFFNESS_MATRIX !<Is .TRUE. if the stiffness matrix is in the equations mapping, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_1",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has already been finished.",err,error,*999)
      ELSE
        EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
          CASE(EQUATIONS_STATIC)
            CALL FlagError("Can not set dynamic matrices for static equations.",err,error,*999)
          CASE(EQUATIONS_QUASISTATIC)
            CALL FlagError("Can not set dynamic matrices for quasi-static equations.",err,error,*999)
          CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
            IF(.NOT.DAMPING_MATRIX) CALL FLAG_WARNING("No damping matrix for first order dynamic equations.",err,error,*999)
            CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_ALL(EQUATIONS_MAPPING,.FALSE.,DAMPING_MATRIX,STIFFNESS_MATRIX, &
              err,error,*999)
          CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
            CALL FlagError("Need to specify three matrices to set for second order dynamic equations.",err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations time dependence type of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",err,error))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_1")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_1",err,error)
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_1")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_1

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrices involved in a second order dynamic equations mapping
  SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_2(EQUATIONS_MAPPING,MASS_MATRIX,DAMPING_MATRIX,STIFFNESS_MATRIX, &
    & err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set the first order matrices for
    LOGICAL, INTENT(IN) :: MASS_MATRIX !<Is .TRUE. if the mass matrix is in the equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: DAMPING_MATRIX !<Is .TRUE. if the damping matrix is in the equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: STIFFNESS_MATRIX !<Is .TRUE. if the stiffness matrix is in the equations mapping, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_2",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has already been finished.",err,error,*999)
      ELSE
        EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
          CASE(EQUATIONS_STATIC)
            CALL FlagError("Can not set dynamic matrices for static equations.",err,error,*999)
          CASE(EQUATIONS_QUASISTATIC)
            CALL FlagError("Can not set dynamic matrices for quasi-static equations.",err,error,*999)
          CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
            IF(MASS_MATRIX) THEN
              CALL FlagError("The mass matrix cannot be present for first order dynamic equations.",err,error,*999)
            ELSE
              IF(.NOT.DAMPING_MATRIX) CALL FLAG_WARNING("No damping matrix for a first order dynamic system.",err,error,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_ALL(EQUATIONS_MAPPING,.FALSE.,DAMPING_MATRIX,STIFFNESS_MATRIX, &
                err,error,*999)
            ENDIF
          CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
            IF(.NOT.MASS_MATRIX) CALL FLAG_WARNING("No mass matrix for a second order dynamic system.",err,error,*999)
            CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_ALL(EQUATIONS_MAPPING,MASS_MATRIX,DAMPING_MATRIX, &
              & STIFFNESS_MATRIX,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations time dependence type of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",err,error))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_2")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_2",err,error)
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_2")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET_2

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrix coefficients in a first order dynamic equations mapping
  SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1(EQUATIONS_MAPPING,DAMPING_MATRIX_COEFFICIENT, &
    & STIFFNESS_MATRIX_COEFFICIENT,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set 
    REAL(DP), INTENT(IN) :: DAMPING_MATRIX_COEFFICIENT !<The damping matrix coefficient
    REAL(DP), INTENT(IN) :: STIFFNESS_MATRIX_COEFFICIENT !<The stiffness matrix coefficient
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has already been finished.",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
            CASE(EQUATIONS_STATIC)
              CALL FlagError("Can not set dynamic matrix coefficients for static equations.",err,error,*999)
            CASE(EQUATIONS_QUASISTATIC)
              CALL FlagError("Can not set dynamic matrix coefficients for quasi-static equations.",err,error,*999)
            CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
              IF(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER)= &
                  & STIFFNESS_MATRIX_COEFFICIENT
              ENDIF
              IF(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER)= &
                  & DAMPING_MATRIX_COEFFICIENT
              ENDIF
            CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
              CALL FlagError("Need to specify three matrix coefficients for second order dynamic equations.", &
                & err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1",err,error)
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrix coefficients in a second order dynamic equations mapping
  SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2(EQUATIONS_MAPPING,MASS_MATRIX_COEFFICIENT, &
    & DAMPING_MATRIX_COEFFICIENT,STIFFNESS_MATRIX_COEFFICIENT,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set 
    REAL(DP), INTENT(IN) :: MASS_MATRIX_COEFFICIENT !<The mass matrix coefficient
    REAL(DP), INTENT(IN) :: DAMPING_MATRIX_COEFFICIENT !<The damping matrix coefficient
    REAL(DP), INTENT(IN) :: STIFFNESS_MATRIX_COEFFICIENT !<The stiffness matrix coefficient
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has already been finished.",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
            CASE(EQUATIONS_STATIC)
              CALL FlagError("Can not set dynamic matrices for static equations.",err,error,*999)
            CASE(EQUATIONS_QUASISTATIC)
              CALL FlagError("Can not set dynamic matrices for quasi-static equations.",err,error,*999)
            CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
              CALL FlagError("Need to specify two matrix coefficients for second order dynamic equations.",err,error,*999)
            CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
              IF(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER)= &
                  & STIFFNESS_MATRIX_COEFFICIENT
              ENDIF
              IF(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER)= &
                  & DAMPING_MATRIX_COEFFICIENT
              ENDIF
              IF(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER)= &
                  & MASS_MATRIX_COEFFICIENT
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The equations time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2",err,error)
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set dynamic matrices
  SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,DYNAMIC_VARIABLE_TYPE,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: DYNAMIC_VARIABLE_TYPE !<The variable type associated with the equations set dynamic matrices.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    LOGICAL :: IS_RESIDUAL_TYPE

    CALL Enters("EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping have been finished.",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          IF(DYNAMIC_VARIABLE_TYPE==0) THEN
            CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE=0
          ELSE
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                IF(EQUATIONS%TIME_DEPENDENCE==EQUATIONS_FIRST_ORDER_DYNAMIC.OR. &
                  EQUATIONS%TIME_DEPENDENCE==EQUATIONS_SECOND_ORDER_DYNAMIC) THEN
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN                 
                    !Check the dynamic variable type is not being by other equations matrices or vectors
                    IF(EQUATIONS%LINEARITY==EQUATIONS_NONLINEAR) THEN
                      IS_RESIDUAL_TYPE=.FALSE.
                      DO matrix_idx=1,CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES
                        IF(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES(matrix_idx)==DYNAMIC_VARIABLE_TYPE) THEN
                          IS_RESIDUAL_TYPE=.TRUE.
                        ENDIF
                      ENDDO
                      IF(IS_RESIDUAL_TYPE.NEQV..TRUE.) THEN
                        LOCAL_ERROR="The specified dynamic variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(DYNAMIC_VARIABLE_TYPE,"*",err,error))// &
                          & " is not the same as any residual variable type."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                    END IF
                    IF(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE==DYNAMIC_VARIABLE_TYPE) THEN
                      LOCAL_ERROR="The specified dynamic variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(DYNAMIC_VARIABLE_TYPE,"*",err,error))// &
                        & " is the same as the variable type for the RHS vector."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                    DO matrix_idx=1,CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                      IF(CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)==DYNAMIC_VARIABLE_TYPE) THEN
                        LOCAL_ERROR="The specified dynamic variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(DYNAMIC_VARIABLE_TYPE,"*",err,error))// &
                          & " is the same as the variable type for linear matrix number "// &
                          & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",err,error))//"."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                    ENDDO !matrix_idx
                    !Check the dynamic variable type is defined on the dependent field
                    IF(DYNAMIC_VARIABLE_TYPE>=1.AND.DYNAMIC_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                      IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(DYNAMIC_VARIABLE_TYPE)%PTR)) THEN
                        EQUATIONS_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE=DYNAMIC_VARIABLE_TYPE
                      ELSE
                        LOCAL_ERROR="The specified dynamic variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(DYNAMIC_VARIABLE_TYPE,"*",err,error))// &
                          & " is not defined on the dependent field."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The specified dynamic variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(DYNAMIC_VARIABLE_TYPE,"*",err,error))// &
                        & " is invalid. The number must either be zero or >= 1 and <= "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Dependent field is not associated",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The equations are not dynamic equations.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations equations set is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations Jacobian map and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE(EQUATIONS_JACOBIAN_TO_VAR_MAP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_JACOBIAN_TO_VAR_MAP_TYPE) :: EQUATIONS_JACOBIAN_TO_VAR_MAP !<The equations Jacobian to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE",err,error,*999)
    
    IF(ALLOCATED(EQUATIONS_JACOBIAN_TO_VAR_MAP%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP)) &
      & DEALLOCATE(EQUATIONS_JACOBIAN_TO_VAR_MAP%EQUATIONS_COLUMN_TO_DOF_VARIABLE_MAP)
    
    CALL Exits("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE",err,error)    
    CALL Exits("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a variable to equations Jacobian map.
  SUBROUTINE EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE(EQUATIONS_JACOBIAN_TO_VAR_MAP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_JACOBIAN_TO_VAR_MAP_TYPE) :: EQUATIONS_JACOBIAN_TO_VAR_MAP !<The equations Jacobian to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE",err,error,*999)
    
    EQUATIONS_JACOBIAN_TO_VAR_MAP%VARIABLE_TYPE=0
    NULLIFY(EQUATIONS_JACOBIAN_TO_VAR_MAP%VARIABLE)
    NULLIFY(EQUATIONS_JACOBIAN_TO_VAR_MAP%JACOBIAN)
    EQUATIONS_JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS=0
    EQUATIONS_JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT=0
    NULLIFY(EQUATIONS_JACOBIAN_TO_VAR_MAP%COLUMN_DOFS_MAPPING)    
    
    CALL Exits("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE",err,error)    
    CALL Exits("EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise an equations matrix to variable maps and deallocate all memory.
  SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE(EQUATIONS_MATRIX_TO_VAR_MAP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TO_VAR_MAP_TYPE) :: EQUATIONS_MATRIX_TO_VAR_MAP !<The equations matrix to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE",err,error,*999)

    IF(ALLOCATED(EQUATIONS_MATRIX_TO_VAR_MAP%COLUMN_TO_DOF_MAP)) &
      & DEALLOCATE(EQUATIONS_MATRIX_TO_VAR_MAP%COLUMN_TO_DOF_MAP)
    
    CALL Exits("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE",err,error)    
    CALL Exits("EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an equations matrix to variable maps.
  SUBROUTINE EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE(EQUATIONS_MATRIX_TO_VAR_MAP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TO_VAR_MAP_TYPE) :: EQUATIONS_MATRIX_TO_VAR_MAP !<The equations matrix to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE",err,error,*999)

    EQUATIONS_MATRIX_TO_VAR_MAP%MATRIX_NUMBER=0
    EQUATIONS_MATRIX_TO_VAR_MAP%VARIABLE_TYPE=0
    NULLIFY(EQUATIONS_MATRIX_TO_VAR_MAP%VARIABLE)
    EQUATIONS_MATRIX_TO_VAR_MAP%NUMBER_OF_COLUMNS=0
    EQUATIONS_MATRIX_TO_VAR_MAP%MATRIX_COEFFICIENT=1.0_DP !Matrices in an equation set are added by default
    NULLIFY(EQUATIONS_MATRIX_TO_VAR_MAP%COLUMN_DOFS_MAPPING)
    
    CALL Exits("EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE",err,error)    
    CALL Exits("EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_EQUATS_MATRIX_TO_VAR_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_FINALISE(EQUATIONS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_FINALISE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
       !Row dofs mappings are linked to the field mapping therefore do not deallocate here
       NULLIFY(EQUATIONS_MAPPING%ROW_DOFS_MAPPING)
       CALL EQUATIONS_MAPPING_DYNAMIC_MAPPING_FINALISE(EQUATIONS_MAPPING%DYNAMIC_MAPPING,err,error,*999)
       CALL EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE(EQUATIONS_MAPPING%LINEAR_MAPPING,err,error,*999)
       CALL EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE(EQUATIONS_MAPPING%NONLINEAR_MAPPING,err,error,*999)
       CALL EquationsMapping_LhsMappingFinalise(EQUATIONS_MAPPING%lhsMapping,err,error,*999)      
       CALL EQUATIONS_MAPPING_RHS_MAPPING_FINALISE(EQUATIONS_MAPPING%RHS_MAPPING,err,error,*999)      
       CALL EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE(EQUATIONS_MAPPING%SOURCE_MAPPING,err,error,*999)      
       CALL EQUATIONS_MAPPING_CREATE_VALUES_CACHE_FINALISE(EQUATIONS_MAPPING%CREATE_VALUES_CACHE,err,error,*999)
       DEALLOCATE(EQUATIONS_MAPPING)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_FINALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_FINALISE")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_INITIALISE(EQUATIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to initialise the equations mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL Enters("EQUATIONS_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(ASSOCIATED(EQUATIONS%EQUATIONS_MAPPING)) THEN
        CALL FlagError("Equations mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS%EQUATIONS_MAPPING,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equations equations mapping.",err,error,*999)
        EQUATIONS%EQUATIONS_MAPPING%EQUATIONS=>EQUATIONS
        EQUATIONS%EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED=.FALSE.
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%ROW_DOFS_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%DYNAMIC_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%LINEAR_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%SOURCE_MAPPING)
        NULLIFY(EQUATIONS%EQUATIONS_MAPPING%CREATE_VALUES_CACHE)
        CALL EQUATIONS_MAPPING_CREATE_VALUES_CACHE_INITIALISE(EQUATIONS%EQUATIONS_MAPPING,err,error,*999)        
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_FINALISE(EQUATIONS%EQUATIONS_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL Errors("EQUATIONS_MAPPING_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_INITIALISE")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set residual vector.
  SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_VARIABLES_NUMBER_SET(EQUATIONS_MAPPING,NUMBER_OF_VARIABLES,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_VARIABLES !<The number of residual variables for this equations set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: PREVIOUS_NUMBER,MIN_NUMBER
    INTEGER(INTG), ALLOCATABLE :: NEW_RESIDUAL_VARIABLE_TYPES(:)
    REAL(DP), ALLOCATABLE :: NEW_RESIDUAL_COEFFICIENTS(:)
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE

    CALL Enters("EQUATIONS_MAPPING_RESIDUAL_VARIABLES_NUMBER_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping have been finished.",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          PREVIOUS_NUMBER=CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES
          IF(NUMBER_OF_VARIABLES/=PREVIOUS_NUMBER) THEN
            MIN_NUMBER=MIN(NUMBER_OF_VARIABLES,PREVIOUS_NUMBER)
            !Create new residual_variable_types array and copy over previous values
            ALLOCATE(NEW_RESIDUAL_VARIABLE_TYPES(NUMBER_OF_VARIABLES),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new residual variable types.",err,error,*999)
            NEW_RESIDUAL_VARIABLE_TYPES=0
            NEW_RESIDUAL_VARIABLE_TYPES(1:MIN_NUMBER)=CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES(1:MIN_NUMBER)
            CALL MOVE_ALLOC(NEW_RESIDUAL_VARIABLE_TYPES,CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES)
            !Create new residual coefficients array and copy over previous values
            ALLOCATE(NEW_RESIDUAL_COEFFICIENTS(NUMBER_OF_VARIABLES),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new residual variable types.",err,error,*999)
            NEW_RESIDUAL_COEFFICIENTS=1.0_DP
            NEW_RESIDUAL_COEFFICIENTS(1:MIN_NUMBER)=CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENTS(1:MIN_NUMBER)
            CALL MOVE_ALLOC(NEW_RESIDUAL_VARIABLE_TYPES,CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES)
            !Set number of residual variables
            CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES=NUMBER_OF_VARIABLES
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF

    CALL Exits("EQUATIONS_MAPPING_RESIDUAL_VARIABLES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(NEW_RESIDUAL_VARIABLE_TYPES)) DEALLOCATE(NEW_RESIDUAL_VARIABLE_TYPES)
    IF(ALLOCATED(NEW_RESIDUAL_COEFFICIENTS)) DEALLOCATE(NEW_RESIDUAL_COEFFICIENTS)
    CALL Errors("EQUATIONS_MAPPING_RESIDUAL_VARIABLES_NUMBER_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_RESIDUAL_VARIABLES_NUMBER_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_VARIABLES_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping LHS mapping and deallocates all memory
  SUBROUTINE EquationsMapping_LhsMappingFinalise(lhsMapping, err, error,*)

    !Argument variables
    TYPE(EquationsMapping_LhsType), POINTER :: lhsMapping !<A pointer to the LHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    CALL Enters("EquationsMapping_LhsMappingFinalise",err,error,*999)

    IF(ASSOCIATED(lhsMapping)) THEN
      IF(ALLOCATED(lhsMapping%LhsDofToEquationsRowMap)) DEALLOCATE(lhsMapping%LhsDofToEquationsRowMap)
      IF(ALLOCATED(lhsMapping%EquationsRowToLhsDofMap)) DEALLOCATE(lhsMapping%EquationsRowToLhsDofMap)
      DEALLOCATE(lhsMapping)
    ENDIF
       
    CALL Exits("EquationsMapping_LhsMappingFinalise")
    RETURN
999 CALL Errors("EquationsMapping_LhsMappingFinalise",err,error)
    CALL Exits("EquationsMapping_LhsMappingFinalise")
    RETURN 1
  END SUBROUTINE EquationsMapping_LhsMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping LHS mapping
  SUBROUTINE EquationsMapping_LhsMappingInitialise(equationsMapping, err, error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping !<A pointer to the equations mapping to initialise the LHS mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    CALL Enters("EquationsMapping_LhsMappingInitialise",err,error,*998)

    IF(ASSOCIATED(equationsMapping)) THEN
      IF(ASSOCIATED(equationsMapping%LhsMapping)) THEN
        CALL FlagError("Equations mapping LHS mapping is already associated.",err,error,*998)
      ELSE
        ALLOCATE(equationsMapping%LhsMapping,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equations mapping RHS mapping.",err,error,*999)
        equationsMapping%LhsMapping%EquationsMapping=>equationsMapping        
        equationsMapping%LhsMapping%LhsVariableType=0
        NULLIFY(equationsMapping%LhsMapping%LhsVariable)
        NULLIFY(equationsMapping%LhsMapping%LhsVariableMapping)
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated.",err,error,*998)
    ENDIF
       
    CALL Exits("EquationsMapping_LhsMappingInitialise")
    RETURN
999 CALL EquationsMapping_LhsMappingFinalise(equationsMapping%LhsMapping,dummyErr,dummyError,*998)
998 CALL Errors("EquationsMapping_LhsMappingInitialise",err,error)
    CALL Exits("EquationsMapping_LhsMappingInitialise")
    RETURN 1
  END SUBROUTINE EquationsMapping_LhsMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set LHS.
  SUBROUTINE EquationsMapping_LhsVariableTypeSet(equationsMapping, lhsVariableType, err, error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: lhsVariableType !<The variable type associated with the equations set LHS. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,vectorIdx
    LOGICAL :: dynamicVariableMatch,linearVariableMatch,nonlinearVariableMatch
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: createValuesCache
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(VARYING_STRING) :: localError

    CALL Enters("EquationsMapping_LhsVariableTypeSet",err,error,*999)

    IF(ASSOCIATED(equationsMapping)) THEN
      IF(equationsMapping%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has been finished.",err,error,*999)
      ELSE
        createValuesCache=>equationsMapping%CREATE_VALUES_CACHE
        IF(ASSOCIATED(createValuesCache)) THEN
          equations=>equationsMapping%Equations
            IF(ASSOCIATED(equations)) THEN
              equationsSet=>equations%EQUATIONS_SET
              IF(ASSOCIATED(equationsSet)) THEN
                dependentField=>equationsSet%Dependent%DEPENDENT_FIELD
                IF(ASSOCIATED(dependentField)) THEN
                  !Check the LHS variable type is being by other equations matrices or vectors
                  dynamicVariableMatch=.FALSE.
                  nonlinearVariableMatch=.FALSE.
                  linearVariableMatch=.FALSE.
                  dynamicVariableMatch=createValuesCache%DYNAMIC_VARIABLE_TYPE==lhsVariableType
                  DO vectorIdx=1,createValuesCache%NUMBER_OF_RESIDUAL_VARIABLES
                    nonlinearVariableMatch=nonlinearVariableMatch.AND. &
                      & createValuesCache%RESIDUAL_VARIABLE_TYPES(vectorIdx)==lhsVariableType
                  ENDDO !vectorIdx
                  DO matrixIdx=1,createValuesCache%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    linearVariableMatch=linearVariableMatch.AND. &
                      & createValuesCache%LINEAR_MATRIX_VARIABLE_TYPES(matrixIdx)==lhsVariableType
                  ENDDO !matrixIdx
                  IF(dynamicVariableMatch.OR.nonlinearVariableMatch.OR.linearVariableMatch) THEN
                    !Check the LHS variable type is defined on the dependent field
                    IF(lhsVariableType>=1.AND.lhsVariableType<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                      IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(lhsVariableType)%Ptr)) THEN
                        createValuesCache%LhsVariableType=lhsVariableType
                      ELSE
                        localError="The specified LHS variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(lhsVariableType,"*",err,error))// &
                          & " is not defined on the dependent field."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                  ELSE
                    localError="The specified LHS variable type of "//TRIM(NUMBER_TO_VSTRING(lhsVariableType,"*",err,error))// &
                      & " is not mapped to any dynamic or linear matrices or nonlinear vectors."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Dependent field is not associated",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations equations set is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF
       
    CALL Exits("EquationsMapping_LhsVariableTypeSet")
    RETURN
999 CALL Errors("EquationsMapping_LhsVariableTypeSet",err,error)
    CALL Exits("EquationsMapping_LhsVariableTypeSet")
    RETURN 1
  END SUBROUTINE EquationsMapping_LhsVariableTypeSet
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping linear mapping and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE(LINEAR_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING !<A pointer to the linear mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,variable_type
 
    CALL Enters("EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE",err,error,*999)

    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
      IF(ALLOCATED(LINEAR_MAPPING%LINEAR_MATRIX_VARIABLE_TYPES)) DEALLOCATE(LINEAR_MAPPING%LINEAR_MATRIX_VARIABLE_TYPES)
      IF(ALLOCATED(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS)) THEN
        DO variable_type=1,SIZE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS,1)
          CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE(LINEAR_MAPPING% &
            & VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type),err,error,*999)
        ENDDO !variable_type
        DEALLOCATE(LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS)        
      ENDIF
      IF(ALLOCATED(LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS)) THEN
        DO matrix_idx=1,SIZE(LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS,1)
          CALL EQUATIONS_MAPPING_EQUATIONS_MATRIX_TO_VAR_MAP_FINALISE(LINEAR_MAPPING% &
            & EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx),err,error,*999)
        ENDDO !matrix_idx
        DEALLOCATE(LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS)
      ENDIF
      DEALLOCATE(LINEAR_MAPPING)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping linear mapping
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE(EQUATIONS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the linear mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL Enters("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%LINEAR_MAPPING)) THEN
        CALL FlagError("Equations mapping linear mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%LINEAR_MAPPING,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equations mapping linear mapping.",err,error,*999)
        EQUATIONS_MAPPING%LINEAR_MAPPING%EQUATIONS_MAPPING=>EQUATIONS_MAPPING       
        EQUATIONS_MAPPING%LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
        EQUATIONS_MAPPING%LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES=0
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_LINEAR_MAPPING_FINALISE(EQUATIONS_MAPPING%LINEAR_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL Errors("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the linear equations matrices in an equation set. 
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFS_SET(EQUATIONS_MAPPING,LINEAR_MATRIX_COEFFICIENTS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping.
    REAL(DP), INTENT(IN) :: LINEAR_MATRIX_COEFFICIENTS(:) !<LINEAR_MATRIX_COEFFICIENTS(matrixIdx). The linear matrix coefficient for the matrixIdx'th linear matrix.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFS_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping is finished.",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN          
          IF(SIZE(LINEAR_MATRIX_COEFFICIENTS,1)==CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES) THEN
            CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS=LINEAR_MATRIX_COEFFICIENTS          
          ELSE
            LOCAL_ERROR="Invalid size of linear matrix coefficeints. The size of the supplied array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(LINEAR_MATRIX_COEFFICIENTS,1),"*",err,error))// &
              & ") must match the number of linear equations matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",err,error))//")."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFS_SET")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFS_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFS_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_COEFFS_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of linear equations matrices
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,NUMBER_OF_LINEAR_EQUATIONS_MATRICES,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set the number of matrices for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_LINEAR_EQUATIONS_MATRICES !<The number of linear equations matrices for the mapping.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    INTEGER(INTG), ALLOCATABLE :: OLD_LINEAR_MATRIX_VARIABLE_TYPES(:)
    REAL(DP), ALLOCATABLE :: OLD_LINEAR_MATRIX_COEFFICIENTS(:)
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has been finished",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
            IF(ASSOCIATED(EQUATIONS_SET)) THEN            
              !Check number of matrices to create is valid
              SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
              CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
                SELECT CASE(EQUATIONS%LINEARITY)
                CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                  IF(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE==0) THEN                  
                    IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES<1) THEN
                      LOCAL_ERROR="The specified number of linear matrices of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",err,error))// &
                        & " is invalid. For non-dynamic linear problems without a equations set RHS the number must be "// &
                        & ">= 1."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ELSE
                    IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES<1) THEN
                      LOCAL_ERROR="The specified number of linear matrices of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",err,error))// &
                        & " is invalid. For non-dynamic linear problems with a equations set RHS the number "// &
                        & "must be >= 1."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ENDIF
                CASE(EQUATIONS_NONLINEAR)
                  IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES<0.OR. &
                    & NUMBER_OF_LINEAR_EQUATIONS_MATRICES>FIELD_NUMBER_OF_VARIABLE_TYPES-2) THEN
                    LOCAL_ERROR="The specified number of linear matrices of "// &
                      & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",err,error))// &
                      & ") is invalid. For non-dynamic non-linear problems the number must be between >= 0 and <= "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES-2,"*",err,error))
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  ENDIF
                CASE DEFAULT
                  LOCAL_ERROR="The equations linearity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",err,error))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
                SELECT CASE(EQUATIONS%LINEARITY)
                CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                  IF(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE==0) THEN                  
                    IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES<1.OR. &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES>FIELD_NUMBER_OF_VARIABLE_TYPES-1) THEN
                      LOCAL_ERROR="The specified number of linear matrices of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",err,error))// &
                        & " is invalid. For dynamic linear problems without a equations set RHS the number must be "// &
                        & "between >= 1 and <= "//TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES-1,"*",err,error))
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ELSE
                    IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES<0) THEN
                      LOCAL_ERROR="The specified number of linear matrices of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",err,error))// &
                        & " is invalid. For dynamic linear problems with a equations set RHS the number "// &
                        & "must be >= 0."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ENDIF
                CASE(EQUATIONS_NONLINEAR)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations linearity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",err,error))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT              
              CASE DEFAULT
                LOCAL_ERROR="The equations time dependence type of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",err,error))//" is invalid."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
              !If we need to reallocate and reset all the create_values cache arrays and change the number of matrices
              IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES/=CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES) THEN
                IF(CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN                  
                  ALLOCATE(OLD_LINEAR_MATRIX_VARIABLE_TYPES(CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate old linear matrix variable types.",err,error,*999)
                  ALLOCATE(OLD_LINEAR_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate old linear matrix coefficients.",err,error,*999)
                  OLD_LINEAR_MATRIX_VARIABLE_TYPES=CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES
                  OLD_LINEAR_MATRIX_COEFFICIENTS=CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS
                ENDIF
                IF(ALLOCATED(CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES)) &
                  & DEALLOCATE(CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES)
                IF(ALLOCATED(CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS)) &
                  & DEALLOCATE(CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS)
                ALLOCATE(CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate linear matrix variable types.",err,error,*999)
                ALLOCATE(CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate linear matrix coefficients.",err,error,*999)
                IF(CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>0) THEN                  
                  IF(NUMBER_OF_LINEAR_EQUATIONS_MATRICES>CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES) THEN
                    CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(1:CREATE_VALUES_CACHE% &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES)=OLD_LINEAR_MATRIX_VARIABLE_TYPES
                    CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(CREATE_VALUES_CACHE% &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES+1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)= &
                      & OLD_LINEAR_MATRIX_VARIABLE_TYPES(1)
                    CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE% &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES)=OLD_LINEAR_MATRIX_COEFFICIENTS
                    CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE% &
                      & NUMBER_OF_LINEAR_EQUATIONS_MATRICES+1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)= &
                      & OLD_LINEAR_MATRIX_COEFFICIENTS(1)
                  ELSE
                    CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)= &
                      & OLD_LINEAR_MATRIX_VARIABLE_TYPES(1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)
                    CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)= &
                      & OLD_LINEAR_MATRIX_COEFFICIENTS(1:NUMBER_OF_LINEAR_EQUATIONS_MATRICES)
                  ENDIF
                ELSE
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                    CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES=0
                    SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
                    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
                      SELECT CASE(EQUATIONS%LINEARITY)
                      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                        IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR)) THEN
                          CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(1)=DEPENDENT_FIELD% &
                            & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%VARIABLE_TYPE
                        ELSE
                          CALL FlagError("Not implemented.",err,error,*999)
                        ENDIF
                        DO matrix_idx=2,CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+1)%PTR)) THEN
                            CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)= &
                              & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+1)%PTR%VARIABLE_TYPE
                          ELSE
                            CALL FlagError("Not implemented.",err,error,*999)
                          ENDIF
                        ENDDO !matrix_idx
                      CASE(EQUATIONS_NONLINEAR)
                        DO matrix_idx=1,CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+2)%PTR)) THEN
                            CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)= &
                              & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+2)%PTR%VARIABLE_TYPE
                          ELSE
                            CALL FlagError("Not implemented.",err,error,*999)
                          ENDIF
                        ENDDO !matrix_idx
                      CASE DEFAULT
                        LOCAL_ERROR="The equations linearity type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",err,error))//" is invalid."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      END SELECT
                    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
                      SELECT CASE(EQUATIONS%LINEARITY)
                      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                        DO matrix_idx=1,CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+2)%PTR)) THEN
                            CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)= &
                              & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+2)%PTR%VARIABLE_TYPE
                          ELSE
                            CALL FlagError("Not implemented.",err,error,*999)
                          ENDIF
                        ENDDO !matrix_idx
                      CASE(EQUATIONS_NONLINEAR)
                        CALL FlagError("Not implemented.",err,error,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The equations linearity type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",err,error))//" is invalid."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      END SELECT
                    CASE DEFAULT
                      LOCAL_ERROR="The equations time dependence type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",err,error))//" is invalid."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    END SELECT
                    CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS=1.0_DP !Equations matrices are added by default
                  ELSE
                    CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                  ENDIF
                ENDIF
                CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                IF(ALLOCATED(OLD_LINEAR_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_LINEAR_MATRIX_VARIABLE_TYPES)
                IF(ALLOCATED(OLD_LINEAR_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_LINEAR_MATRIX_COEFFICIENTS)
              ENDIF
            ELSE
              CALL FlagError("Equations equations set is not associated",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations mapping equations is not associated",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_LINEAR_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_LINEAR_MATRIX_VARIABLE_TYPES)    
    IF(ALLOCATED(OLD_LINEAR_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_LINEAR_MATRIX_COEFFICIENTS)    
    CALL Errors("EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Sets the mapping between the dependent field variable types and the linear equations matrices
  SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,LINEAR_MATRIX_VARIABLE_TYPES,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping
    INTEGER(INTG), INTENT(IN) :: LINEAR_MATRIX_VARIABLE_TYPES(:) !<The matrix variable types to map to each linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL Enters("EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has been finished.",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          IF(SIZE(LINEAR_MATRIX_VARIABLE_TYPES,1)==CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES) THEN
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  !Check input values
                  DO matrix_idx=1,CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    IF(LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)/=0) THEN
                      !Check the residual variable type is not being by other equations matrices or vectors
                      !Don't check against the residual variable as we can have linear parts of nonlinear equations
                      IF(CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE==LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)) THEN
                        LOCAL_ERROR="The specified linear matrix variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx),"*",err,error))// &
                          & " for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",err,error))// &
                          & " is the same as the variable type for the dynamic matrices."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                      IF(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE==LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)) THEN
                        LOCAL_ERROR="The specified linear matrix variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx),"*",err,error))// &
                          & " for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",err,error))// &
                          & " is the same as the variable type for the RHS vector."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF                       
                      !Check to see if the linear matrix variable numbers are defined on the dependent field
                      IF(LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)>=1.OR. &
                        & LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                        IF(.NOT.ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx))%PTR)) THEN
                          LOCAL_ERROR="The linear matrix variable type of "// &
                            & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx),"*",err,error))// &
                            & " for linear matrix NUMBER "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",err,error))// &
                            & " is not defined on the dependent field."
                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="The linear matrix variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx),"*",err,error))// &
                          & " for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",err,error))// &
                          & " is invalid. The variable types must be either zero or >= 1 and <= "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                    ENDIF
                  ENDDO !matrix_idx
                  CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES=LINEAR_MATRIX_VARIABLE_TYPES
                ELSE
                  CALL FlagError("Dependent field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations equations set is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="Invalid size of linear matrix variable types. The size of the supplied array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(LINEAR_MATRIX_VARIABLE_TYPES,1),"*",err,error))// &
              & ") must match the number of linear equations matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,"*",err,error))//")."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated.",err,error,*999)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping nonlinear mapping and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE(NONLINEAR_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING !<A pointer to the nonlinear mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) matrix_idx
 
    CALL Enters("EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE",err,error,*999)

    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
      DO matrix_idx=1,NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
        CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE(NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP(matrix_idx),err,error,*999)
        CALL EQUATIONS_MAPPING_EQUATS_JACOBIAN_TO_VAR_MAP_FINALISE(NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(matrix_idx),err,error,*999)
      ENDDO
      IF(ALLOCATED(NONLINEAR_MAPPING%EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP)) &
        & DEALLOCATE(NONLINEAR_MAPPING%EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP)
      IF(ALLOCATED(NONLINEAR_MAPPING%RESIDUAL_VARIABLES)) &
        & DEALLOCATE(NONLINEAR_MAPPING%RESIDUAL_VARIABLES)
      IF(ALLOCATED(NONLINEAR_MAPPING%RESIDUAL_COEFFICIENTS)) &
        & DEALLOCATE(NONLINEAR_MAPPING%RESIDUAL_COEFFICIENTS)
      IF(ALLOCATED(NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP)) &
        & DEALLOCATE(NONLINEAR_MAPPING%VAR_TO_JACOBIAN_MAP)
      IF(ALLOCATED(NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP)) &
        & DEALLOCATE(NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP)
      DEALLOCATE(NONLINEAR_MAPPING)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping nonlinear mapping
  SUBROUTINE EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE(EQUATIONS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the nonlinear mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL Enters("EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%NONLINEAR_MAPPING)) THEN
        CALL FlagError("Equations mapping nonlinear mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%NONLINEAR_MAPPING,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equations mapping nonlinear mapping.",err,error,*999)
        EQUATIONS_MAPPING%NONLINEAR_MAPPING%EQUATIONS_MAPPING=>EQUATIONS_MAPPING
        EQUATIONS_MAPPING%NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES=0
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL Exits("EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_NONLINEAR_MAPPING_FINALISE(EQUATIONS_MAPPING%NONLINEAR_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL Errors("EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_NONLINEAR_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set residual vector.
  SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_COEFFS_SET(EQUATIONS_MAPPING,RESIDUAL_COEFFICIENTS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: RESIDUAL_COEFFICIENTS(:) !<RESIDUAL_COEFFICIENTS(residualIdx). The coefficient applied to the residualIdx'th residual vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_RESIDUAL_COEFFS_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping have been finished.",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          IF(SIZE(RESIDUAL_COEFFICIENTS,1)==CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES) THEN
            CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENTS=RESIDUAL_COEFFICIENTS          
          ELSE
            LOCAL_ERROR="Invalid size of residual coefficients. The size of the supplied array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(RESIDUAL_COEFFICIENTS,1),"*",err,error))// &
              & ") must match the number of residuals ("// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MAPPING%CREATE_VALUES_CACHE% &
              & NUMBER_OF_RESIDUAL_VARIABLES,"*",err,error))//")."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_RESIDUAL_COEFFS_SET")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_RESIDUAL_COEFFS_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_RESIDUAL_COEFFS_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_COEFFS_SET

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set residual vector.
  SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,RESIDUAL_VARIABLE_TYPES,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: RESIDUAL_VARIABLE_TYPES(:) !<RESIDUAL_VARIABLE_TYPE(variable_idx). The variable_idx'th variable type associated with the equations set residual vector. The first variable type must correspond to the diagonal terms in the full solver Jacobian so that the solver mapping can use boundary conditions on this first variable to decide whether to keep rows.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,variable_idx,NUMBER_OF_RESIDUAL_VARIABLES,RESIDUAL_VARIABLE_TYPE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping have been finished.",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          NUMBER_OF_RESIDUAL_VARIABLES=SIZE(RESIDUAL_VARIABLE_TYPES,1)
          IF(NUMBER_OF_RESIDUAL_VARIABLES==CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES) THEN
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                IF(EQUATIONS%LINEARITY==EQUATIONS_NONLINEAR) THEN
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                    !Check the residual variable types are not being used by other equations matrices or vectors
                    DO variable_idx=1,NUMBER_OF_RESIDUAL_VARIABLES
                      RESIDUAL_VARIABLE_TYPE=RESIDUAL_VARIABLE_TYPES(variable_idx)
                      IF(EQUATIONS%TIME_DEPENDENCE==EQUATIONS_STATIC .OR. EQUATIONS%TIME_DEPENDENCE==EQUATIONS_QUASISTATIC) THEN
                        IF(CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE==RESIDUAL_VARIABLE_TYPE) THEN
                          LOCAL_ERROR="The specified residual variable type of "// &
                            & TRIM(NUMBER_TO_VSTRING(RESIDUAL_VARIABLE_TYPE,"*",err,error))// &
                            & " is the same as the variable type for the dynamic matrices."
                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                        ENDIF
                      ELSE IF(EQUATIONS%TIME_DEPENDENCE==EQUATIONS_FIRST_ORDER_DYNAMIC.OR. & 
                        & EQUATIONS%TIME_DEPENDENCE==EQUATIONS_SECOND_ORDER_DYNAMIC) THEN
                        IF(CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE/=RESIDUAL_VARIABLE_TYPE) THEN
                          LOCAL_ERROR="The specified residual variable type of "// &
                            & TRIM(NUMBER_TO_VSTRING(RESIDUAL_VARIABLE_TYPE,"*",err,error))// &
                            & " is not the same as the variable type for the dynamic matrices."
                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("The equations set time dependence is not set.",err,error,*999)
                      END IF
                      IF(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE==RESIDUAL_VARIABLE_TYPE) THEN
                        LOCAL_ERROR="The specified residual variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(RESIDUAL_VARIABLE_TYPE,"*",err,error))// &
                          & " is the same as the variable type for the RHS vector."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                      DO matrix_idx=1,CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                        IF(CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)==RESIDUAL_VARIABLE_TYPE) THEN
                          LOCAL_ERROR="The specified residual variable type of "// &
                            & TRIM(NUMBER_TO_VSTRING(RESIDUAL_VARIABLE_TYPE,"*",err,error))// &
                            & " is the same as the variable type for linear matrix number "// &
                            & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",err,error))//"."
                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                        ENDIF
                      ENDDO !matrix_idx
                      !Check the residual variable number is defined on the dependent field
                      IF(RESIDUAL_VARIABLE_TYPE>=1.AND.RESIDUAL_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                        IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(RESIDUAL_VARIABLE_TYPE)%PTR)) THEN
                          CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES(variable_idx)=RESIDUAL_VARIABLE_TYPE
                        ELSE
                          LOCAL_ERROR="The specified residual variable type of "// &
                            & TRIM(NUMBER_TO_VSTRING(RESIDUAL_VARIABLE_TYPE,"*",err,error))// &
                            & " is not defined on the dependent field."
                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="The specified residual variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(RESIDUAL_VARIABLE_TYPE,"*",err,error))// &
                          & " is invalid. The variable type must either be zero or >= 1 and <= "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                    ENDDO !variable_idx
                  ELSE
                    CALL FlagError("Dependent field is not associated",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The equations set is not a nonlinear equations set.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations equations set is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="Invalid number of variables. The number of residual variables " &
              & //TRIM(NUMBER_TO_VSTRING(NUMBER_OF_RESIDUAL_VARIABLES,"*",err,error)) &
              & //" should be "//TRIM(NUMBER_TO_VSTRING(CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES,"*",err,error))
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET
  
  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set RHS vector.
  SUBROUTINE EQUATIONS_MAPPING_RHS_COEFF_SET(EQUATIONS_MAPPING,RHS_COEFFICIENT,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: RHS_COEFFICIENT!<The coefficient applied to the equations set RHS vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_RHS_COEFF_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE/=0) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT=RHS_COEFFICIENT
          ELSE
            CALL FlagError("The equations mapping RHS variable type has not been set.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_RHS_COEFF_SET")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_RHS_COEFF_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_RHS_COEFF_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RHS_COEFF_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping RHS mapping and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_RHS_MAPPING_FINALISE(RHS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING !<A pointer to the RHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL Enters("EQUATIONS_MAPPING_RHS_MAPPING_FINALISE",err,error,*999)

    IF(ASSOCIATED(RHS_MAPPING)) THEN
      IF(ALLOCATED(RHS_MAPPING%RHS_DOF_TO_EQUATIONS_ROW_MAP)) DEALLOCATE(RHS_MAPPING%RHS_DOF_TO_EQUATIONS_ROW_MAP)
      IF(ALLOCATED(RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP)) DEALLOCATE(RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP)
      DEALLOCATE(RHS_MAPPING)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_RHS_MAPPING_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_RHS_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RHS_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping RHS mapping
  SUBROUTINE EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE(EQUATIONS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the RHS mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL Enters("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%RHS_MAPPING)) THEN
        CALL FlagError("Equations mapping RHS mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%RHS_MAPPING,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equations mapping RHS mapping.",err,error,*999)
        EQUATIONS_MAPPING%RHS_MAPPING%EQUATIONS_MAPPING=>EQUATIONS_MAPPING        
        EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE_TYPE=0
        NULLIFY(EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE)
        NULLIFY(EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE_MAPPING)
        EQUATIONS_MAPPING%RHS_MAPPING%RHS_COEFFICIENT=1.0_DP
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_RHS_MAPPING_FINALISE(EQUATIONS_MAPPING%RHS_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL Errors("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RHS_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set rhs vector.
  SUBROUTINE EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,RHS_VARIABLE_TYPE,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: RHS_VARIABLE_TYPE !<The variable type associated with the equations set rhs vector. If the problem does not have a rhs vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has been finished.",err,error,*999)
      ELSE
        CREATE_VALUES_CACHE=>EQUATIONS_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          IF(RHS_VARIABLE_TYPE==0) THEN
            CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=0
          ELSE
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  !Check the RHS variable type is not being by other equations matrices or vectors
                  IF(CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE==RHS_VARIABLE_TYPE) THEN
                    LOCAL_ERROR="The specified RHS variable type of "// &
                      & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",err,error))// &
                        & " is the same as the variable type for the dynamic matrices."
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  ENDIF
                  DO matrix_idx=1,CREATE_VALUES_CACHE%NUMBER_OF_RESIDUAL_VARIABLES
                    IF(CREATE_VALUES_CACHE%RESIDUAL_VARIABLE_TYPES(matrix_idx)==RHS_VARIABLE_TYPE) THEN
                      LOCAL_ERROR="The specified RHS variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",err,error))// &
                        & " is the same as the variable type for the residual vector."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ENDDO
                  DO matrix_idx=1,CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    IF(CREATE_VALUES_CACHE%LINEAR_MATRIX_VARIABLE_TYPES(matrix_idx)==RHS_VARIABLE_TYPE) THEN
                      LOCAL_ERROR="The specified RHS variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",err,error))// &
                        & " is the same as the variable type for linear matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",err,error))//"."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ENDDO !matrix_idx
                  !Check the RHS variable number is defined on the dependent field
                  IF(RHS_VARIABLE_TYPE>=1.AND.RHS_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                    IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(RHS_VARIABLE_TYPE)%PTR)) THEN
                      CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE=RHS_VARIABLE_TYPE
                    ELSE
                      LOCAL_ERROR="The specified RHS variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",err,error))// &
                        & " is not defined on the dependent field."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The specified RHS variable type of "//TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",err,error))// &
                      & " is invalid. The number must either be zero or >= 1 and <= "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Dependent field is not associated",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations equations set is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set source vector.
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_COEFF_SET(EQUATIONS_MAPPING,SOURCE_COEFFICIENT,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: SOURCE_COEFFICIENT!<The coefficient applied to the equations set source vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_SOURCE_COEFF_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping has been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE/=0) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_COEFFICIENT=SOURCE_COEFFICIENT
          ELSE
            CALL FlagError("The equations mapping source variable type has not been set.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_SOURCE_COEFF_SET")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_SOURCE_COEFF_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_SOURCE_COEFF_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_COEFF_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping source mapping and deallocates all memory
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE(SOURCE_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING !<A pointer to the SOURCE mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL Enters("EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE",err,error,*999)

    IF(ASSOCIATED(SOURCE_MAPPING)) THEN
      DEALLOCATE(SOURCE_MAPPING)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping source mapping
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE(EQUATIONS_MAPPING,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to initialise the source mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL Enters("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(ASSOCIATED(EQUATIONS_MAPPING%SOURCE_MAPPING)) THEN
        CALL FlagError("Equations mapping source mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_MAPPING%SOURCE_MAPPING,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equations mapping source mapping.",err,error,*999)
        EQUATIONS_MAPPING%SOURCE_MAPPING%EQUATIONS_MAPPING=>EQUATIONS_MAPPING        
        EQUATIONS_MAPPING%SOURCE_MAPPING%SOURCE_VARIABLE_TYPE=0
        NULLIFY(EQUATIONS_MAPPING%SOURCE_MAPPING%SOURCE_VARIABLE)
        EQUATIONS_MAPPING%SOURCE_MAPPING%SOURCE_COEFFICIENT=1.0_DP
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL Exits("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE")
    RETURN
999 CALL EQUATIONS_MAPPING_SOURCE_MAPPING_FINALISE(EQUATIONS_MAPPING%SOURCE_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL Errors("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE",err,error)
    CALL Exits("EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a source field variable and the equations set source vector.
  SUBROUTINE EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,SOURCE_VARIABLE_TYPE,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: SOURCE_VARIABLE_TYPE !<The variable type associated with the equations set source vector. If the problem does not have a source vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL Enters("EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        CALL FlagError("Equations mapping have been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(SOURCE_VARIABLE_TYPE==0) THEN
            EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE=0
          ELSE
            EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
                  SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
                  IF(ASSOCIATED(SOURCE_FIELD)) THEN                    
                    !Check the source variable type is defined on the source field
                    IF(SOURCE_VARIABLE_TYPE>=1.AND.SOURCE_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                      IF(ASSOCIATED(SOURCE_FIELD%VARIABLE_TYPE_MAP(SOURCE_VARIABLE_TYPE)%PTR)) THEN
                        EQUATIONS_MAPPING%CREATE_VALUES_CACHE%SOURCE_VARIABLE_TYPE=SOURCE_VARIABLE_TYPE
                      ELSE
                        LOCAL_ERROR="The specified source variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(SOURCE_VARIABLE_TYPE,"*",err,error))// &
                          & " is not defined on the source field."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The specified source variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(SOURCE_VARIABLE_TYPE,"*",err,error))// &
                        & " is invalid. The number must either be zero or >= 1 and <= "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Source field is not associated",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set source is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations equations set is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations mapping equations is not associated.",err,error,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FlagError("Equations mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations mapping is not associated",err,error,*999)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET",err,error)
    CALL Exits("EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Finalise an equations mapping equations matrix map.
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE(VAR_TO_EQUATIONS_COLUMN_MAP,err,error,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_COLUMN_MAP_TYPE) :: VAR_TO_EQUATIONS_COLUMN_MAP !<The variable dof to equations column map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE",err,error,*999)
    
    IF(ALLOCATED(VAR_TO_EQUATIONS_COLUMN_MAP%COLUMN_DOF)) &
      & DEALLOCATE(VAR_TO_EQUATIONS_COLUMN_MAP%COLUMN_DOF)
    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE",err,error)    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations Jacobian map and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE(VAR_TO_EQUATIONS_JACOBIAN_MAP,err,error,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_JACOBIAN_MAP_TYPE) :: VAR_TO_EQUATIONS_JACOBIAN_MAP !<The variable to equations Jacobian map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE",err,error,*999)
    
    IF(ALLOCATED(VAR_TO_EQUATIONS_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP)) &
      & DEALLOCATE(VAR_TO_EQUATIONS_JACOBIAN_MAP%DOF_TO_COLUMNS_MAP)
    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE",err,error)    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a variable to equations Jacobian map
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE(VAR_TO_EQUATIONS_JACOBIAN_MAP,err,error,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_JACOBIAN_MAP_TYPE) :: VAR_TO_EQUATIONS_JACOBIAN_MAP !<The variable to equations Jacobian map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE",err,error,*999)
    
    VAR_TO_EQUATIONS_JACOBIAN_MAP%VARIABLE_TYPE=0
    NULLIFY(VAR_TO_EQUATIONS_JACOBIAN_MAP%VARIABLE)
    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE",err,error)    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_JACOBIAN_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations matrices map and deallocates all memory.
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE(VAR_TO_EQUATIONS_MATRICES_MAP,err,error,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_MATRICES_MAP_TYPE) :: VAR_TO_EQUATIONS_MATRICES_MAP !<The variable to equations matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL Enters("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE",err,error,*999)
    
    IF(ALLOCATED(VAR_TO_EQUATIONS_MATRICES_MAP%EQUATIONS_MATRIX_NUMBERS)) &
      & DEALLOCATE(VAR_TO_EQUATIONS_MATRICES_MAP%EQUATIONS_MATRIX_NUMBERS)
    IF(ALLOCATED(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS)) THEN
      DO matrix_idx=1,SIZE(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS,1)
        CALL EQUATIONS_MAPPING_VAR_TO_EQUATS_COLUMN_MAP_FINALISE(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS( &
          & matrix_idx),err,error,*999)
      ENDDO !matrix_idx
      DEALLOCATE(VAR_TO_EQUATIONS_MATRICES_MAP%DOF_TO_COLUMNS_MAPS)
    ENDIF
    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE",err,error)    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an equations mapping equations matrix map.
  SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE(VAR_TO_EQUATIONS_MATRICES_MAP,err,error,*)

    !Argument variables
    TYPE(VAR_TO_EQUATIONS_MATRICES_MAP_TYPE) :: VAR_TO_EQUATIONS_MATRICES_MAP !<The variable to equations matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE",err,error,*999)

    VAR_TO_EQUATIONS_MATRICES_MAP%VARIABLE_INDEX=0
    VAR_TO_EQUATIONS_MATRICES_MAP%VARIABLE_TYPE=0
    NULLIFY(VAR_TO_EQUATIONS_MATRICES_MAP%VARIABLE)
    VAR_TO_EQUATIONS_MATRICES_MAP%NUMBER_OF_EQUATIONS_MATRICES=0
    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE")
    RETURN
999 CALL Errors("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE",err,error)    
    CALL Exits("EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE

  !
  !================================================================================================================================
  !
  
END MODULE EQUATIONS_MAPPING_ROUTINES

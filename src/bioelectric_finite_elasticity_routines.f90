!> \file
!> \authors Thomas Heidlauf
!> \brief This module handles all routines pertaining to bioelectrics coupled with finite elasticity.
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

!>This module handles all routines pertaining to bioelectrics coupled with finite elasticity.


MODULE BIOELECTRIC_FINITE_ELASTICITY_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BIOELECTRIC_ROUTINES
  USE BIODOMAIN_EQUATION_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_IO_ROUTINES
  USE FIELD_ROUTINES
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TYPES

  IMPLICIT NONE

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET
  
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_POST_LOOP
  
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_UPDATE_GEOMETRIC_FIELD

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a bioelectrics finite elasticity equation type of a multi physics equations set class.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectrics finite elasticity equation type of a multi physics equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity equation.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string


    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    CALL FLAG_ERROR("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP is not implemented.",ERR,ERROR,*999)

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a bioelectrics finite elasticity equation finite element equations set.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    CALL FLAG_ERROR("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE is not implemented.",ERR,ERROR,*999)

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a bioelectrics finite elasticity equation type of a multi physics equations set class.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)

    CALL FLAG_ERROR("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET is not implemented.",ERR,ERROR,*999)

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a bioelectric finite elasticity problem type .
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE
      CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectric finite elasticity problem type of a multi physics problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectric finite elasticity problem.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: MONODOMAIN_SUB_LOOP,ELASTICITY_SUB_LOOP
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS,MONODOMAIN_SOLVERS,ELASTICITY_SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(MONODOMAIN_SUB_LOOP)
    NULLIFY(ELASTICITY_SUB_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVERS)
    NULLIFY(MONODOMAIN_SOLVERS)
    NULLIFY(ELASTICITY_SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(CELLML_EQUATIONS)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)

      !--------------------------------------------------------------------
      !   Transient Gudunov monodomain, simple finite elasticity  
      !--------------------------------------------------------------------
      CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(CONTROL_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up the control sub loop for monodomain
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(MONODOMAIN_SUB_LOOP,'MONODOMAIN_TIME_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(MONODOMAIN_SUB_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(MONODOMAIN_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up the control sub loop for finite elasicity
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(ELASTICITY_SUB_LOOP,'ELASTICITY_LOAD_INCREMENT_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(ELASTICITY_SUB_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(ELASTICITY_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
            !Sub-loops are finished when parent is finished
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the monodomain sub loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(MONODOMAIN_SOLVERS,2,ERR,ERROR,*999)
            !Set the first solver to be a differential-algebraic equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"ODE Solver",ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_RESTART_SET(SOLVER,.TRUE.,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)

            !Get the finite elasticity sub loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(ELASTICITY_SOLVERS,1,ERR,ERROR,*999)
            !Set the finite elasticity solver to be a nonlinear solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the monodomain solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(MONODOMAIN_SOLVERS,ERR,ERROR,*999)

            !Get the finite elasticity solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(ELASTICITY_SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a bioelectrics finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)

            !Get the monodomain sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Create the solver equations for the second (parabolic) solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)

            !Get the finite elasticity sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Get the finite elasticity solver and create the finite elasticity solver equations
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            
            !Get the monodomain sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Get the solver equations for the second (parabolic) solver
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             

            !Get the finite elasticity sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Finish the creation of the finite elasticity solver equations
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Create the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Get the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectrics finite elasticity equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a transient monodomain quasistatic finite elasticity equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the bioelectrics finite elasticity problem pre-solve.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              CALL BIODOMAIN_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
              CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem post solve.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)
            SELECT CASE(SOLVER%SOLVE_TYPE)
            CASE(SOLVER_DAE_TYPE)
              CALL BIOELECTRIC_POST_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_TYPE)
              CALL BIOELECTRIC_POST_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_NONLINEAR_TYPE)
              CALL FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Solver solve type "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration???
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          SELECT CASE(PROBLEM%TYPE)
          CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              !do nothing ???
            CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
              IF(CONTROL_LOOP%PROBLEM%SUBTYPE==PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE) THEN
                CALL BIOELECTRIC_FIN_ELA_INDEPENDENT_FIELD_INTERPOLATE(CONTROL_LOOP,ERR,ERROR,*999)
              ENDIF
              CALL FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for bioelectric finite elasticity problem type."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
              & " is not valid for a multi physics problem class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - do nothing!
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG) :: equations_set_idx
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(REGION_TYPE), POINTER :: DEPENDENT_REGION   
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: FILENAME,LOCAL_ERROR,METHOD
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ELASTICITY_SUB_LOOP,BIOELECTRIC_SUB_LOOP

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            SELECT CASE(PROBLEM%TYPE)
            CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
              !the monodomain time loop - output of the monodomain fields
              CALL BIODOMAIN_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
                & " is not valid for a multi physics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
            CALL BIOELECTRIC_FINITE_ELASTICITY_UPDATE_GEOMETRIC_FIELD(CONTROL_LOOP,.FALSE.,ERR,ERROR,*999)
          CASE DEFAULT
            !do nothing
          END SELECT
        ELSE
          CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - output the finite elasticity fields 
        IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          !Export the dependent field for this time step
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            PROBLEM=>CONTROL_LOOP%PROBLEM
            IF(ASSOCIATED(PROBLEM)) THEN
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(ELASTICITY_SUB_LOOP)
              !Get the solver. The first solver of the second sub loop will contain the finite elasticity dependent field equation set
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              !Loop over the equations sets associated with the solver
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  DO equations_set_idx=1,SOLVER_MAPPING%numberOfEquationsSets
                    EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(equations_set_idx)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      NULLIFY(DEPENDENT_REGION)
                      CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                      FILENAME="MainTime_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                        & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                      METHOD="FORTRAN"
                      CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="Equations set is not associated for equations set index "// &
                        & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                        & " in the solver mapping."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                ELSE
                  CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver solver equations are not associated.",ERR,ERROR,*999)
              ENDIF
              IF(PROBLEM%SUBTYPE==PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE) THEN
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(BIOELECTRIC_SUB_LOOP)
                !Get the solver. The second solver of the first sub loop will contain the bioelectrics equation set
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,BIOELECTRIC_SUB_LOOP,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(BIOELECTRIC_SUB_LOOP,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                !Loop over the equations sets associated with the solver
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    DO equations_set_idx=1,SOLVER_MAPPING%numberOfEquationsSets
                      EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(equations_set_idx)%ptr
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        NULLIFY(DEPENDENT_REGION)
                        CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                        FILENAME="MainTime_M_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                          & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                        METHOD="FORTRAN"
                        CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                        
                        WRITE(*,*) TIME_LOOP%ITERATION_NUMBER
                        
                      ELSE
                        LOCAL_ERROR="Equations set is not associated for equations set index "// &
                          & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                          & " in the solver mapping."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !equations_set_idx
                  ELSE
                    CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver solver equations are not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF !PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE
            ELSE
              CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_POST_LOOP")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_POST_LOOP",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_POST_LOOP")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

  !>Update the the bioelectric equation geometric field from the finite elasticity dependent field (deformed geometry)
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_UPDATE_GEOMETRIC_FIELD(CONTROL_LOOP,CALC_CLOSEST_GAUSS_POINT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    LOGICAL, INTENT(IN) :: CALC_CLOSEST_GAUSS_POINT !<If true then the closest finite elasticity Gauss point for each bioelectrics node is calculated
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT,CONTROL_LOOP_ELASTICITY,CONTROL_LOOP_MONODOMAIN
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_ELASTICITY,GEOMETRIC_FIELD_MONODOMAIN,GEOMETRIC_FIELD_ELASTICITY
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_MONODOMAIN,INDEPENDENT_FIELD_MONODOMAIN,DEPENDENT_FIELD_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING,NODES_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,FIELD_VARIABLE_M,FIELD_VARIABLE_FE
    INTEGER(INTG) :: component_idx,element_idx,ne
    INTEGER(INTG) :: DEPENDENT_FIELD_INTERPOLATION,GEOMETRIC_FIELD_INTERPOLATION
    INTEGER(INTG) :: node_idx,node_idx_2,NODE_LEFT,NODE_RIGHT,NUMBER_OF_NODES,GAUSS_POINT,gauss_idx,fibre_idx
    INTEGER(INTG) :: nodes_in_Xi_1,nodes_in_Xi_2,nodes_in_Xi_3,n3,n2,n1,dof_idx,dof_idx2
    REAL(DP) :: XVALUE_M,XVALUE_FE,DIST_LEFT,DIST_RIGHT,VALUE,VALUE_LEFT,VALUE_RIGHT,DISTANCE
    REAL(DP) :: XI(3),PREVIOUS_NODE(3)
    LOGICAL :: OUTSIDE_NODE
    REAL(DP), POINTER :: GAUSS_POSITIONS(:,:)

    CALL ENTERS("BIOELECTRIC_FINITE_ELASTICITY_UPDATE_GEOMETRIC_FIELD",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP_ROOT)
    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_ELASTICITY)
    NULLIFY(CONTROL_LOOP_MONODOMAIN)
    NULLIFY(PROBLEM)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(INDEPENDENT_FIELD_ELASTICITY)
    NULLIFY(DEPENDENT_FIELD_MONODOMAIN)
    NULLIFY(INDEPENDENT_FIELD_MONODOMAIN)
    NULLIFY(DEPENDENT_FIELD_ELASTICITY)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(GEOMETRIC_FIELD_MONODOMAIN)
    NULLIFY(GEOMETRIC_FIELD_ELASTICITY)
    NULLIFY(ELEMENTS_MAPPING)
    NULLIFY(INTERPOLATED_POINT)
    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(FIELD_VARIABLE_M)
    NULLIFY(FIELD_VARIABLE_FE)
    NULLIFY(GAUSS_POSITIONS)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          SELECT CASE(PROBLEM%TYPE)
          CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
            SELECT CASE(PROBLEM%SUBTYPE)

            CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)

              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
              !get the monodomain sub loop, solvers, solver, and finally geometric and field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    GEOMETRIC_FIELD_MONODOMAIN=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_MONODOMAIN)) THEN
                      CALL FLAG_ERROR("Geometric field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_MAPPING)
              NULLIFY(EQUATIONS_SET)
              NULLIFY(SOLVER_EQUATIONS)
              !get the finite elasticity sub loop, solvers, solver, and finally the dependent field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    DEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(DEPENDENT_FIELD_ELASTICITY)) THEN
                      CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              DO component_idx=1,GEOMETRIC_FIELD_MONODOMAIN%VARIABLES(1)%NUMBER_OF_COMPONENTS
                !check for identical interpolation of the fields
                GEOMETRIC_FIELD_INTERPOLATION=GEOMETRIC_FIELD_MONODOMAIN%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
                DEPENDENT_FIELD_INTERPOLATION=DEPENDENT_FIELD_ELASTICITY%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
                IF(GEOMETRIC_FIELD_INTERPOLATION==DEPENDENT_FIELD_INTERPOLATION) THEN
                  !copy the dependent field components to the geometric field
                  CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & component_idx,ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="The interpolation type of component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR, &
                    & ERROR))//" of field number "//TRIM(NUMBER_TO_VSTRING(GEOMETRIC_FIELD_MONODOMAIN%USER_NUMBER,"*",ERR, &
                    & ERROR))//" does not coincide with the interpolation type of field number " &
                    & //TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD_ELASTICITY%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO

            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)

              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
              !get the monodomain sub loop, solvers, solver, and finally geometric field and dependent field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    GEOMETRIC_FIELD_MONODOMAIN=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_MONODOMAIN)) THEN
                      CALL FLAG_ERROR("Geometric field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    ! the Field_V_Variable_Type contains the 3D nodal positions
                    DEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(DEPENDENT_FIELD_MONODOMAIN)) THEN
                      CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_MONODOMAIN)) THEN
                      CALL FLAG_ERROR("Independent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_MAPPING)
              NULLIFY(EQUATIONS_SET)
              NULLIFY(SOLVER_EQUATIONS)
              !get the finite elasticity sub loop, solvers, solver, and finally the dependent and independent fields
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    INDEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_ELASTICITY)) THEN
                      CALL FLAG_ERROR("Independent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    GEOMETRIC_FIELD_ELASTICITY=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_ELASTICITY)) THEN
                      CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              
              node_idx=0
              fibre_idx=0
              CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE_M,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VARIABLE_FE,ERR,ERROR,*999)

              ELEMENTS_MAPPING=>GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION% &
                & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
              NODES_MAPPING=>GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION% &
                & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES
              
              !loop through the elements of the finite elasticity mesh (internal and boundary elements)
              !no need to consider ghost elements here since only bioelectrical fields are changed
              DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%BOUNDARY_FINISH
                ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)

                !get the finite elasticity dependent field interpolation parameters of this element
                INTERPOLATION_PARAMETERS=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS &
                  & (FIELD_U_VARIABLE_TYPE)%PTR
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
                INTERPOLATED_POINT=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR

                !the Field_V_Variable_Type of the FE independent field contains the number of nodes in each Xi-direction of the bioelectrics grid
                dof_idx=FIELD_VARIABLE_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,nodes_in_Xi_1,ERR,ERROR,*999)
                dof_idx=FIELD_VARIABLE_FE%COMPONENTS(2)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,nodes_in_Xi_2,ERR,ERROR,*999)
                dof_idx=FIELD_VARIABLE_FE%COMPONENTS(3)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,nodes_in_Xi_3,ERR,ERROR,*999)

                !if there is no bioelectrics grid in this finite elasticity element, jump to the next finite elasticity element
                IF((nodes_in_Xi_1==0).OR.(nodes_in_Xi_2==0).OR.(nodes_in_Xi_3==0)) CYCLE
                
                !get the positions of the Gauss points of the Finite Elasticity element, GAUSS_POSITIONS(components,number_of_Gauss_points)
                GAUSS_POSITIONS=>GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION% &
                  & MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP( &
                  & BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_POSITIONS
                
                !assume Xi(1) to be normal to the seed surface, i.e. the seed points have Xi(1)=0
                XI=[0.0_DP,1.0_DP/(REAL(2*nodes_in_Xi_2)),1.0_DP/(REAL(2*nodes_in_Xi_3))]
                
                !assume that the bioelectrics node numbers are increased in order Xi(1), Xi(2), Xi(3) 
                DO n3=1,nodes_in_Xi_3
                  DO n2=1,nodes_in_Xi_2
                    fibre_idx=fibre_idx+1
                    DO n1=1,nodes_in_Xi_1
                      node_idx=node_idx+1
                      
                      !store the fibre number this bioelectrics node belongs to.
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1,node_idx,3,fibre_idx,ERR,ERROR,*999) 

                      !find the interpolated position of the bioelectric grid node from the FE dependent field
                      CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
                      !update the bioelectrics dependent field Field_V_Variable_Type
                      !the Field_V_Variable_Type of the monodomain dependent field contains the nodal positions in 3D
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
                        & INTERPOLATED_POINT%VALUES(1,1),ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1,node_idx,2, &
                        & INTERPOLATED_POINT%VALUES(2,1),ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1,node_idx,3, &
                        & INTERPOLATED_POINT%VALUES(3,1),ERR,ERROR,*999)

                      IF(n1==1) THEN
                        !a new line of bioelectrics grid nodes begins
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,0.0_DP,ERR,ERROR,*999)
                      ELSE
                        !get the position in 3D of the previous node
                        dof_idx2=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                          & DERIVATIVES(1)%VERSIONS(1)
                        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,dof_idx2,PREVIOUS_NODE(1),ERR,ERROR,*999)
                        dof_idx2=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                          & DERIVATIVES(1)%VERSIONS(1)
                        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,dof_idx2,PREVIOUS_NODE(2),ERR,ERROR,*999)
                        dof_idx2=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                          & DERIVATIVES(1)%VERSIONS(1)
                        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,dof_idx2,PREVIOUS_NODE(3),ERR,ERROR,*999)

                        !compute the distance between the previous node and the actual node
                        VALUE=SQRT( &
                          & (INTERPOLATED_POINT%VALUES(1,1)-PREVIOUS_NODE(1))*(INTERPOLATED_POINT%VALUES(1,1)-PREVIOUS_NODE(1))+ &
                          & (INTERPOLATED_POINT%VALUES(2,1)-PREVIOUS_NODE(2))*(INTERPOLATED_POINT%VALUES(2,1)-PREVIOUS_NODE(2))+ &
                          & (INTERPOLATED_POINT%VALUES(3,1)-PREVIOUS_NODE(3))*(INTERPOLATED_POINT%VALUES(3,1)-PREVIOUS_NODE(3)))

                        !get the position in 1D of the previous node
                        dof_idx2=FIELD_VARIABLE_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx-1)% &
                          & DERIVATIVES(1)%VERSIONS(1)
                        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,dof_idx2,VALUE_LEFT,ERR,ERROR,*999)
                        !update the current 1D node position
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,VALUE_LEFT+VALUE,ERR,ERROR,*999)
                      ENDIF
                        
                      IF(CALC_CLOSEST_GAUSS_POINT) THEN
                        !calculate the closest finite elasticity Gauss point of each bioelectrics node
                        DISTANCE=1000000.0_DP
                        GAUSS_POINT=0
                        DO gauss_idx=1,SIZE(GAUSS_POSITIONS,2)
                          !compute the distance between the bioelectrics node and the Gauss point
                          VALUE=SQRT( &
                            & (Xi(1)-GAUSS_POSITIONS(1,gauss_idx))*(Xi(1)-GAUSS_POSITIONS(1,gauss_idx))+ &
                            & (Xi(2)-GAUSS_POSITIONS(2,gauss_idx))*(Xi(2)-GAUSS_POSITIONS(2,gauss_idx))+ &
                            & (Xi(3)-GAUSS_POSITIONS(3,gauss_idx))*(Xi(3)-GAUSS_POSITIONS(3,gauss_idx)))
                          IF(VALUE<DISTANCE) THEN
                            DISTANCE=VALUE
                            GAUSS_POINT=gauss_idx
                          ENDIF
                        ENDDO !gauss_idx
                        IF(GAUSS_POINT==0) CALL FLAG_WARNING("Closest Gauss Point not found",ERR,ERROR,*999)
                        !store the nearest Gauss Point info and the inElement info (local element number!!!)
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,4,GAUSS_POINT,ERR,ERROR,*999)
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,1,node_idx,5,ne,ERR,ERROR,*999)
                      ENDIF !CALC_CLOSEST_GAUSS_POINT
                      
                      XI(1)=XI(1)+1.0_DP/(REAL(nodes_in_Xi_1-1))
                    ENDDO !n1
                    XI(1)=0.0_DP
                    XI(2)=XI(2)+1.0_DP/(REAL(nodes_in_Xi_2))
                  ENDDO !n2
                  XI(1)=0.0_DP
                  XI(2)=1.0_DP/(REAL(2*nodes_in_Xi_2))
                  XI(3)=Xi(3)+1.0_DP/(REAL(nodes_in_Xi_3))
                ENDDO !n3

              ENDDO !element_idx

            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity problem type of a multi physics problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
              & " is not valid for a multi physics problem class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - do nothing!
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_UPDATE_GEOMETRIC_FIELD")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FINITE_ELASTICITY_UPDATE_GEOMETRIC_FIELD",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FINITE_ELASTICITY_UPDATE_GEOMETRIC_FIELD")
    RETURN 1
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_UPDATE_GEOMETRIC_FIELD

  !
  !================================================================================================================================
  !

  !>Interpolates the finite elasticity independent field from the biolectrics independent field.
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available
  SUBROUTINE BIOELECTRIC_FIN_ELA_INDEPENDENT_FIELD_INTERPOLATE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT,CONTROL_LOOP_ELASTICITY,CONTROL_LOOP_MONODOMAIN
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_MONODOMAIN,INDEPENDENT_FIELD_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING,NODES_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE_U,FIELD_VARIABLE_V,FIELD_VARIABLE_FE
    INTEGER(INTG) :: INDEPENDENT_FIELD_INTERPOLATION,GEOMETRIC_FIELD_INTERPOLATION
    INTEGER(INTG) :: node_idx,node_idx_2,NODE,equation_set_idx,element_idx,gauss_idx,ne
    INTEGER(INTG) :: nearestGP,inElement,gp_count,dof_idx
    INTEGER(INTG) :: NUMBER_OF_GAUSS_POINTS
    REAL(DP) :: DISTANCE,XVALUE_FE,XVALUE_M,VALUE
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_GAUSS_POINTS=64
    INTEGER(INTG) :: NUMBER_OF_NODES(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: ACTIVE_STRESS_VALUES(MAX_NUMBER_OF_GAUSS_POINTS)

    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_MONODOMAIN)
    NULLIFY(CONTROL_LOOP_ELASTICITY)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(FIELD_VARIABLE_U)
    NULLIFY(FIELD_VARIABLE_V)
    NULLIFY(FIELD_VARIABLE_FE)
    
    CALL ENTERS("BIOELECTRIC_FIN_ELA_INDEPENDENT_FIELD_INTERPOLATE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      PROBLEM=>CONTROL_LOOP%PROBLEM
      IF(ASSOCIATED(PROBLEM)) THEN
        SELECT CASE(PROBLEM%SUBTYPE)
        CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)
          IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
            !--- MONODOMAIN ---
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_MONODOMAIN)) CALL FLAG_ERROR("Independent field is not associated.", &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF

            !--- FINITE ELASTICITY ---
            NULLIFY(SOLVERS)
            NULLIFY(SOLVER)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_ELASTICITY)) CALL FLAG_ERROR("Independent field is not associated.",ERR, &
                    & ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF

            !--- NOW INTERPOLATE ---
            ELEMENTS_MAPPING=>INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION% &
              & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
            NODES_MAPPING=>INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION% &
              & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES

            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE_U,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VARIABLE_V,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE_FE,ERR,ERROR,*999)

            !loop through the finite elasticity elements
            !first process the internal and boundary elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%BOUNDARY_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              
              NUMBER_OF_GAUSS_POINTS=INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_ELASTICITY% &
                & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP &
                & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS

              IF(NUMBER_OF_GAUSS_POINTS>MAX_NUMBER_OF_GAUSS_POINTS) CALL FLAG_ERROR( & 
                & "NUMBER_OF_GAUSS_POINTS is greater than MAX_NUMBER_OF_GAUSS_POINTS.",ERR,ERROR,*999)
              NUMBER_OF_NODES=0
              ACTIVE_STRESS_VALUES=0.0_DP
              
              !loop through the bioelectrics nodes
              DO node_idx=1,NODES_MAPPING%NUMBER_OF_LOCAL
                dof_idx=FIELD_VARIABLE_V%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,inElement,ERR,ERROR,*999) !component 5 of variable V contains inElem info (LOCAL NUMBERING!!!)

                !check if the bioelectrics node is located within the finite elasticity element
                IF(inElement==ne) THEN
                  !component 4 of variable V contains Nearest Gauss Point info
                  dof_idx=FIELD_VARIABLE_V%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                    & VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & dof_idx,nearestGP,ERR,ERROR,*999)
                  IF(nearestGP>MAX_NUMBER_OF_GAUSS_POINTS) CALL FLAG_ERROR( &
                    & "Nearest Gauss Point is greater than MAX_NUMBER_OF_GAUSS_POINTS.",ERR,ERROR,*999)
                  !component 1 of variable U contains the active stress
                  dof_idx=FIELD_VARIABLE_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                    & VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & dof_idx,VALUE,ERR,ERROR,*999)

                  !count the number of bioelectrics nodes that are closest to each finite elasticity Gauss point
                  NUMBER_OF_NODES(nearestGP)=NUMBER_OF_NODES(nearestGP)+1
                  !add up the active stress value
                  ACTIVE_STRESS_VALUES(nearestGP)=ACTIVE_STRESS_VALUES(nearestGP)+VALUE
                ENDIF
              ENDDO

              !loop throught the finite elasticity Gauss points
              DO gauss_idx=1,NUMBER_OF_GAUSS_POINTS
                !make sure we don't divide by zero
                IF(NUMBER_OF_NODES(gauss_idx)<=0) THEN
                  VALUE=0.0_DP
                ELSE
                  VALUE=ACTIVE_STRESS_VALUES(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                ENDIF
                dof_idx=FIELD_VARIABLE_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
              ENDDO !gauss_idx

            ENDDO !element_idx

            !now the ghost elements -- get the relevant info from the other computational nodes
            CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD_ELASTICITY, & 
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD_ELASTICITY, & 
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="Independent field interpolation is not implemented for problem subtype " &
            & //TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("BIOELECTRIC_FIN_ELA_INDEPENDENT_FIELD_INTERPOLATE")
    RETURN
999 CALL ERRORS("BIOELECTRIC_FIN_ELA_INDEPENDENT_FIELD_INTERPOLATE",ERR,ERROR)
    CALL EXITS("BIOELECTRIC_FIN_ELA_INDEPENDENT_FIELD_INTERPOLATE")
    RETURN 1

  END SUBROUTINE BIOELECTRIC_FIN_ELA_INDEPENDENT_FIELD_INTERPOLATE

  !
  !================================================================================================================================
  !

END MODULE BIOELECTRIC_FINITE_ELASTICITY_ROUTINES

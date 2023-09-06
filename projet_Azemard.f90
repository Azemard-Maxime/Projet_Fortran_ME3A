! Projet Azémard ME3.2
!______________________________________________________________________________________________________!

MODULE Precisions
    IMPLICIT NONE
    INTEGER,PARAMETER :: DP=KIND(1.0d0)
END MODULE Precisions
   





! Exercice 1

PROGRAM projetODC
    USE Precisions
    IMPLICIT NONE

    REAL(KIND = DP) :: eps
    REAL(KIND = DP) :: gamma, r 
    REAL(KIND = DP) :: PB, PH, TB , TH , wi , M_ini_newton , cB , cH , delta_P , Mi_max , Mi_min
    INTEGER :: nbValuesGraph

    ! initialisation des variables de l'exercice

    eps = 1.d-9
    gamma = 1.4d0
    r = 287.d0
    PB = 97825.d0
    TB = 283.1d0
    TH = 298.d0
    wi = 485.d0

    ! recupéré dans le code de l'exercice 1 mais ne peuvent pas être récupérés dans le code de l'exercice 2 si pas d'initialisation
    cB = 0.d0
    cH = 0.d0

    PH = 3200000.d0  ! Attention : PH ici est référencé à la Q2 et est différent de PH dans la Q1
    M_ini_newton = 1.5d0

    delta_P = 350000.d0
    Mi_max = 0.d0

    Mi_min = 1.01d0
    nbValuesGraph = 1000





    !Exercice 1
    WRITE(*,*) " "
    WRITE(*,*) "Exercice 1"
    CALL Exercice1(gamma , r , TB , TH , wi , PB , cB , cH)
    !WRITE(*,*) cB , cH  ! pour vérifier la mise à jour des valeurs de cB et cH

    !Exercice 2
    WRITE(*,*) " "
    WRITE(*,*) "Exercice 2"
    CALL Exercice2(eps, M_ini_newton , gamma, cB , cH , PB ,PH)

    !Exercice 3
    WRITE(*,*) " "
    WRITE(*,*) "Exercice 3"
    CALL Exercice3(eps , PB , gamma , cB , cH , delta_P , Mi_max)

    ! Exercice 4
        ! plot "courbeSurpressionReflechie.dat" , "courbeSurpressionIncidente.dat" , "courbePH.dat" w lp
        ! taper la ligne ci-dessus dans gnuplot pour afficher les courbes
    WRITE(*,*) " "
    WRITE(*,*) "Exercice 4"
    CALL Exercice4(Mi_min , Mi_max , nbValuesGraph, gamma , cB , cH , PB)
    WRITE(*,*) "plot 'courbeSurpressionReflechie.dat' , 'courbeSurpressionIncidente.dat' , 'courbePH.dat' w lp"

    
    ! Exercice 5
    WRITE(*,*) " "
    WRITE(*,*) "Exercice 5"
    CALL Exercice5(eps , PB , gamma , cB , cH)

    ! Exercice 6
    WRITE(*,*) " "
    WRITE(*,*) "Exercice 6"
    CALL Exercice6 (gamma , eps , PB)

END PROGRAM projetODC




!______________________________________________________________________________________________________!



Subroutine Exercice1 (gamma , r , TB , TH , vitesse , PB , cB , cH)
    USE Precisions
    IMPLICIT NONE

    REAL(KIND = DP) :: gamma , r , TB , TH , vitesse , PB , Pi           ! les variables d'entrée
    REAL (KIND = DP) :: celerite , mach , cB , cH , PrH , Pr , surpression_incident , surpression_reflechie                 ! les variables de sortie
    REAL (KIND = DP) , EXTERNAL :: celerite_pression , celerite_temperature , nbMach , pression_Haute , pression_incidente ! les fonctions externes
    REAL (KIND = DP), EXTERNAL :: pression_choc_refl_air , surpression_air ! les fonctions externes

    ! On cherche le nombre de mach
    celerite = celerite_temperature(gamma , r , TB)
    mach = nbMach(vitesse , celerite)
    WRITE(*,*) " M =" , mach

    ! calcul de la pression haute
    cB = celerite_temperature( gamma , r , TB)
    cH = celerite_temperature( gamma , r , TH)
    PrH = pression_Haute ( gamma , mach , cB , cH , PB)
    WRITE(*,*) " PH =" , PrH

    ! calcul de la pression Pi
    Pi = pression_incidente (gamma , mach , PB)
    WRITE(*,*) " Pi =" , Pi

    ! calcul de la pression Pr (PB = P0)
    Pr = pression_choc_refl_air(Pi , PB)
    WRITE(*,*) " Pr =" , Pr

    ! calcul de la surpression
    surpression_incident = surpression_air(Pi , PB)
    WRITE(*,*) " surpression_incident =" , surpression_incident
    surpression_reflechie = surpression_air(Pr , PB)
    WRITE(*,*) " surpression_reflechie =" , surpression_reflechie

END Subroutine Exercice1








Subroutine Exercice2(epsilon , x0 , gamma , CB , CH , PB , PH)
    ! Il faut utiliser newton raphson qui permet de trouver rapidement les racines d'une fonction
    ! la relation PH / PB = ... se met facilement sous la forme d'une fonction f(x) = 0
    ! la seule variable est M le mach du choc
    ! on peut donc utiliser newton raphson pour trouver M

    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: epsilon , x0 , gamma , CB , CH , PB , PH
    REAL(KIND = DP) :: Mi , Pi_Q2
    REAL(KIND = DP) , EXTERNAL :: newton_raphson_Q2 , pression_incidente

    Mi = newton_raphson_Q2(epsilon , x0 , gamma , CB , CH , PB , PH)
    WRITE(*,*) " Mi =" , Mi

    ! calcul de la pression incidente
    Pi_Q2 = pression_incidente (gamma , Mi , PB)
    WRITE(*,*) " Pi_Q2 =" , Pi_Q2

END Subroutine Exercice2









Subroutine Exercice3(epsilon , PB , gamma , CB , CH , delta_P , Mi_max)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: epsilon , PB , gamma , CB , CH , delta_P , Mi_max
    REAL(KIND = DP) :: PI_Max , Wi_max , PH_max , Pr_max , Pi_maxPart2 , x0_Pi , Wi_maxPart2 , mach_max_part2
    REAL (KIND = DP) :: PH_max_Part2
    REAL(KIND = DP) , EXTERNAL :: Vit_choc_inc , nbMach , pression_Haute , newton_raphson_Pi


    ! calcul de la pression incidente maximale
    ! on cherche le mach du choc incident maximal
    ! delta P = 350 kPa


    ! dans cette 1er partie, Pi est calculé rapidement avec delta P et est connu
    PI_max = delta_P + PB
    Wi_max = Vit_choc_inc(CB , gamma , PI_max , PB)
    Mi_max = nbMach(Wi_max , CB)
    WRITE(*,*) " mach_max =" , Mi_max

    ! calcul de PH pour ce mach
    PH_max = pression_Haute ( gamma , Mi_max , CB , CH , PB)
    WRITE(*,*) " PH_max =" , PH_max



    ! dans cette 2eme partie, Pi est calculé par newton raphson
    Pr_max = delta_P + PB
    x0_Pi = 100000.d0
    Pi_maxPart2 = newton_raphson_Pi(epsilon , x0_Pi , PB , Pr_max)

    Wi_maxPart2 = Vit_choc_inc(CB , gamma , Pi_maxPart2 , PB)
    mach_max_part2 = nbMach(Wi_maxPart2 , CB)
    WRITE(*,*) " Mr_max =" , mach_max_part2

    ! calcul de PH pour ce mach
    PH_max_part2 = pression_Haute ( gamma , mach_max_part2 , CB , CH , PB)
    WRITE(*,*) " PHr_max =" , PH_max_part2

END Subroutine Exercice3








Subroutine Exercice4(Mi_min , Mi_max , nbPoints, gamma , CB , CH , PB)
    ! code pour tracer les pressions PH en fonction de Mi entre 1.01 et Mi_max

    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: Mi_min , Mi_max , gamma , CB , CH , PB , pas
    REAL(KIND = DP) , EXTERNAL :: pression_Haute , pression_incidente , surpression_air , pression_choc_refl_air

    REAL(KIND=DP), Allocatable :: monTableau_PH(:,:)
    integer :: i  , nbPoints

    Allocate(monTableau_PH(nbPoints+1,4))
    pas = (Mi_max-Mi_min)/REAL(nbPoints,DP)

    DO i = 1 , nbPoints+1
        monTableau_PH(i,1) = Mi_min + (i-1)*pas
        monTableau_PH(i,2) = pression_Haute ( gamma , monTableau_PH(i,1) , CB , CH , PB)
        monTableau_PH(i,3) = surpression_air (pression_incidente (gamma , monTableau_PH(i,1) , PB) , PB)
        monTableau_PH(i,4) = surpression_air (pression_choc_refl_air (pression_incidente &
        (gamma , monTableau_PH(i,1) , PB) , PB) , PB)
    END DO


    ! on écrit les résultats dans un fichier texte
    OPEN(UNIT=10,FILE="resultats.txt",STATUS="UNKNOWN")
    DO i = 1 , nbPoints+1
        WRITE(10,*) monTableau_PH(i,1) , monTableau_PH(i,2) , monTableau_PH(i,3) , monTableau_PH(i,4)
    END DO
    CLOSE(UNIT=10)

    ! on écrit dans des fichiers dat pour lire les courbes
    OPEN(UNIT=11,FILE="courbePH.dat",STATUS="UNKNOWN")
    DO i = 1 , nbPoints+1
        WRITE(11,*) monTableau_PH(i,1) , monTableau_PH(i,2)
    END DO
    CLOSE(UNIT=11)

    OPEN(UNIT=12,FILE="courbeSurpressionIncidente.dat",STATUS="UNKNOWN")
    DO i = 1 , nbPoints+1
        WRITE(12,*) monTableau_PH(i,1) , monTableau_PH(i,3)
    END DO
    CLOSE(UNIT=12)

    OPEN(UNIT=13,FILE="courbeSurpressionReflechie.dat",STATUS="UNKNOWN")
    DO i = 1 , nbPoints+1
        WRITE(13,*) monTableau_PH(i,1) , monTableau_PH(i,4)
    END DO
    CLOSE(UNIT=13)

    DeAllocate(monTableau_PH)


END Subroutine Exercice4







Subroutine Exercice5(epsilon , PB , gamma , CB , CH)

    USE Precisions
    IMPLICIT NONE

    !REAL (KIND = DP) :: 
    Character (len=100) :: nomFichier
    INTEGER :: i , unit
    REAL (KIND = DP) :: CB , PB , gamma , CH , Pi_ci , Pr_cr , x0_pi_cr , PI_cr , epsilon 
    REAL (KIND = DP) , EXTERNAL :: Vit_choc_inc , pression_Haute , Pi_choc_ref , nbMach , newton_raphson_Pi
    REAL (KIND = DP) , DIMENSION(9,5) :: TabSeuils 

    ! Récupération des valeurs dans le fichiers "seuils.inp"

    nomFichier = "seuils.inp"
    OPEN(UNIT=14,FILE= nomFichier,STATUS="OLD")
    DO i = 1 , 9
        READ(14,*) TabSeuils(i,1) 
        WRITE(*,*) "Pi = " , TabSeuils(i,1)

        ! calcul de la pression incidente pour un choc incident (cf: Ex3)
        Pi_ci = TabSeuils(i,1) + PB
        TabSeuils(i , 3) = nbMach (Vit_choc_inc(CB , gamma , Pi_ci , PB) , CB) ! la vitesse en mach
        WRITE(*,*) "Mi = " , TabSeuils(i,3)
        TabSeuils(i , 2) = pression_Haute ( gamma , TabSeuils(i,3) , CB , CH , PB)  ! calcul de la pression haute
        WRITE(*,*) "PHi = " , TabSeuils(i,2)
        
        ! calcul de la pression incidente pour un choc réfléchi (cf: Ex3)
        Pr_cr = TabSeuils(i,1) + PB
        x0_pi_cr = 100000.d0
        PI_cr = newton_raphson_Pi(epsilon , x0_pi_cr , PB , Pr_cr)
        TabSeuils(i , 5) = nbMach (Vit_choc_inc(CB , gamma , PI_cr , PB) , CB) ! la vitesse en mach
        WRITE(*,*) "Mr = " , TabSeuils(i,5)
        TabSeuils(i , 4) = pression_Haute ( gamma , TabSeuils(i,5) , CB , CH , PB)  ! calcul de la pression haute
        WRITE(*,*) "PHr = " , TabSeuils(i,4)

    END DO
    CLOSE(UNIT=14)

    ! on écrit les résultats dans un fichier "seuil.out"
    unit = 15
    OPEN(UNIT=unit,FILE="seuil.out",STATUS="UNKNOWN")
    WRITE(unit,*) '     Valeur du seuil','   |    ','PH seuil (incident)','   |  ','    Mach Choc       ',&
    '   |   ','  PH seuil (réfléchi)','   |   ','    Mach'
    DO i = 1 , 9
        WRITE(15,*) TabSeuils(i,1) , TabSeuils(i,2) , TabSeuils(i,3) , TabSeuils(i,4) , TabSeuils(i,5)
    END DO

END Subroutine Exercice5



Subroutine Exercice6 (gamma , epsilon , PB)
    USE Precisions
    IMPLICIT NONE

    INTEGER :: choice , i , max_danger_1 , max_danger_2 , max_danger_3 , max_danger_4
    REAL (KIND = DP) :: nbDeMach , ch_pression_Haute , Pi , PB , gamma , cB , cH , epsilon , delta_Pi , delta_Pr , Pr
    REAL (KIND = DP) :: x0_pi , Mi
    REAL (KIND = DP) , EXTERNAL :: pression_incidente , pression_Haute , surpression_air ,&
    newton_raphson_Pi , pression_choc_refl_air , newton_raphson_Q2
    REAL (KIND = DP) , DIMENSION(9,1) :: niveau
    Character(len=100) , DIMENSION (9,1):: textDanger

    
    WRITE (*,*) "choisir nb de mach (1) ou pression_Haute (2) ?"
    READ (*,*) choice

    IF ((choice /= 1) .AND. (choice /= 2)) THEN 
        WRITE(*,*) "erreur de saisie"
        STOP
    END IF

    ! remplir le tableau avec le fichier "danger.inp"
    OPEN(UNIT=16,FILE="danger.inp",STATUS="OLD")
    DO i = 1 , 9
        READ(16,*) niveau(i,1) , textDanger(i,1)
    END DO
    CLOSE(UNIT=16)



    IF (choice == 1) THEN
        WRITE(*,*) "choisir un nombre de mach"
        READ(*,*) nbDeMach

        ! calcul de Pi - PB
        Pi = pression_incidente(gamma , nbDeMach , PB)
        delta_Pi = surpression_air(Pi , PB)


        ! calcul de delta Pr
        Pr = pression_choc_refl_air(Pi , PB)
        delta_Pr = surpression_air(Pr , PB)

        ! comparaison avec les seuils
        WRITE (*,*) " Pour Pi :"
        DO i = 1 , 9
            IF (((delta_Pi - niveau(i,1)) > epsilon) .or. (dabs(delta_Pi - niveau(i,1)) < epsilon))  THEN
                max_Danger_1 = i
            END IF
        END DO
        WRITE(*,*) "danger = " , textDanger(max_Danger_1,1)

        WRITE (*,*) " Pour Pr :"
        DO i = 1 , 9
            IF (((delta_Pr - niveau(i,1)) > epsilon) .or. (dabs(delta_Pr - niveau(i,1)) < epsilon))  THEN
                max_Danger_2 = i
            END IF
        END DO
        WRITE(*,*) "danger = " , textDanger(max_Danger_2,1)

        
    

    ELSE  
        WRITE(*,*) "choisir une pression_Haute"
        READ(*,*) ch_pression_Haute

        ! calcul de Mi
        x0_pi = 100000.d0
        Mi = newton_raphson_Q2(epsilon , x0_pi , gamma , CB , CH , PB , ch_pression_Haute)
        Pi = pression_incidente(gamma , Mi , PB)

        ! calcul de delta Pi
        delta_Pi = surpression_air(Pi , PB)

        ! calcul de delta Pr
        Pr = pression_choc_refl_air(Pi , PB)
        delta_Pr = surpression_air(Pr , PB)


        ! comparaison avec les seuils
        WRITE (*,*) " Pour Pi :"
        DO i = 1 , 9
            IF (((delta_Pi - niveau(i,1)) > epsilon) .or. (dabs(delta_Pi - niveau(i,1)) < epsilon))  THEN
                max_Danger_3 = i
            END IF
        END DO
        !WRITE(*,*) "danger = " ,  textDanger(max_Danger_3,1)
        WRITE(*,*) "Erreur ici : pas trouvé comment faire pour afficher le danger"

        WRITE (*,*) " Pour Pr :"
        DO i = 1 , 9
            IF (((delta_Pr - niveau(i,1)) > epsilon) .or. (dabs(delta_Pr - niveau(i,1)) < epsilon))  THEN
                max_Danger_4 = i
            END IF
        END DO
        !WRITE(*,*) "danger = " , textDanger(max_Danger_4,1)
        WRITE(*,*) "Erreur ici : pas trouvé comment faire pour afficher le danger"

        ! Erreur ici : pas trouvé comment faire pour afficher le danger
    
    END IF


        
END Subroutine Exercice6





!______________________________________________________________________________________________________!
! Fonctions

FUNCTION nbMach (vitesse,celerite)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: vitesse, nbMach , celerite
    nbMach = vitesse / celerite
    RETURN
END FUNCTION nbMach




FUNCTION celerite_pression (gamma,pression,masse_volumique)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: gamma , pression, masse_volumique, celerite_pression
    celerite_pression = DSQRT(gamma * pression / masse_volumique)
END FUNCTION celerite_pression





FUNCTION celerite_temperature (gamma, r , temperature)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: gamma , r , temperature, celerite_temperature
    celerite_temperature = SQRT(gamma * r * temperature)
    RETURN
END FUNCTION celerite_temperature






FUNCTION pression_Haute ( gamma , mach , Cel_B , Cel_H , PB)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: gamma , mach , Cel_B , Cel_H , PB , pression_Haute , partGauche , partDroite , partDroitePuissance
    partGauche = ((gamma - 1.d0) / (gamma + 1.d0))*( (((2.d0*gamma)/(gamma - 1.d0)) * mach**2) - 1.d0)
    partDroite = (1.d0 - (((gamma - 1.d0)/(gamma +1.d0))*(mach**(-1.d0))*(((Cel_B/Cel_H)*(mach**2.d0)) - 1.d0)))
    partDroitePuissance = partDroite ** (- ((2.d0 * gamma) / (gamma - 1.d0)))
    pression_Haute = PB * partGauche * partDroitePuissance
    RETURN
END FUNCTION pression_Haute






FUNCTION pression_incidente (gamma , mach , Pression_B)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: gamma , mach , Pression_B , pression_incidente , rap_gamma , var_A
    rap_gamma = (gamma - 1.d0) / (gamma + 1.d0)
    var_A = ((2.d0 * gamma) / (gamma + 1.d0)) * (mach**2.d0)
    pression_incidente = ( Pression_B * (var_A - rap_gamma))
    RETURN
END FUNCTION pression_incidente







FUNCTION pression_choc_refl_air (Pi , PB)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: Pi , PB , pression_choc_refl_air
    pression_choc_refl_air = (((8.d0 * Pi) - PB) / (Pi + (6.d0 * PB))) * Pi
    RETURN
END FUNCTION pression_choc_refl_air






FUNCTION surpression_air ( PressionChoisie , PB)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: PressionChoisie , PB , surpression_air
    surpression_air = pressionChoisie - PB
    RETURN
END FUNCTION surpression_air







! la fonction suivante est la fonction f(x) = 0
FUNCTION fonction_newton_Mi(gamma , mach , Cel_B , Cel_H , PB , PH)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: gamma , mach , Cel_B , Cel_H , PB , PH 
    REAL(KIND = DP) :: fonction_newton_Mi, var_A , var_B , var_C , rap_gamma

    rap_gamma = (gamma - 1.d0) / (gamma + 1.d0)
    var_A = rap_gamma * ( (((2.d0*gamma)/(gamma - 1.d0)) * mach**2) - 1.d0)
    var_B = (1.d0 - ((rap_gamma)*(mach**(-1.d0))*(((Cel_B/Cel_H)*(mach**2.d0)) - 1.d0)))
    var_C = var_B ** (- ((2.d0 * gamma) / (gamma - 1.d0)))
    fonction_newton_Mi = (var_A * var_C) - (PH / PB)
    RETURN
END FUNCTION fonction_newton_Mi






! la fonction suivante est la dérivée de la fonction f(x) = 0
FUNCTION deriveFonction_newton_Mi(gamma , mach , Cel_B , Cel_H)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: gamma , mach , Cel_B , Cel_H
    REAL(KIND = DP) :: deriveFonction_newton_Mi , rap_gamma

    rap_gamma = (gamma - 1.d0) / (gamma + 1.d0)

    deriveFonction_newton_Mi = rap_gamma*((4.d0*gamma)/(gamma-1.d0))*mach*(1.d0-rap_gamma*(1.d0/mach)&
    *((Cel_B/Cel_H)*(mach**2.d0)-1.d0))**((-2.d0*gamma)/(gamma-1.d0))+rap_gamma*(((2.d0*gamma)/(gamma-1.d0))&
    *(mach**2.d0)-1.d0)*((-2.d0*gamma)/(gamma-1.d0))*(-rap_gamma*(1.d0+1.d0/(mach**2.d0)))&
    *(1.d0-rap_gamma*(1.d0/mach)*((Cel_B/Cel_H)*(mach**2.d0)-1.d0))**((-3.d0*gamma+1.d0)/(gamma-1.d0))
    RETURN
END FUNCTION deriveFonction_newton_Mi






FUNCTION newton_raphson_Q2(epsilon , x0 , gamma , Cel_B , Cel_H , PB , PH)
    USE Precisions
    IMPLICIT NONE
    REAL (KIND = DP) :: epsilon , x0 , gamma , Cel_B , Cel_H , PB , PH
    REAL (KIND = DP) :: newton_raphson_Q2 , x1 , val
    REAL (KIND = DP) , EXTERNAL :: fonction_newton_Mi , deriveFonction_newton_Mi
    Integer :: i

    val = fonction_newton_Mi(gamma , x0 , Cel_B , Cel_H , PB , PH)

    DO WHILE ((DABS(val) > epsilon) .OR. (i < 200))  ! on s'arrête si on a atteint la précision ou si on a fait 200 itérations (pas de convergence)
        val = fonction_newton_Mi(gamma , x0 , Cel_B , Cel_H , PB , PH)
        x1 = x0 - (fonction_newton_Mi(gamma , x0 , Cel_B , Cel_H , PB , PH) / deriveFonction_newton_Mi(gamma , x0 , Cel_B , Cel_H))
        x0 = x1
        i = i + 1
    END DO

    newton_raphson_Q2 = x1
    RETURN
END FUNCTION newton_raphson_Q2






FUNCTION Vit_choc_inc(Cel_B , gamma , Pi , PB)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: Cel_B , gamma , Pi , PB 
    REAL (KIND = DP) :: Vit_choc_inc
    Vit_choc_inc = Cel_B * (SQRT(((gamma - 1.d0)/(2.d0 * gamma)) + (((gamma + 1.d0)/(2.d0 * gamma)) * (Pi / PB))))
    RETURN
END FUNCTION Vit_choc_inc







FUNCTION fonction_newton_Pi(Pression_inc , Pression_B , Pression_refl)
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: Pression_inc , Pression_B , Pression_refl
    REAL(KIND = DP) :: fonction_newton_Pi
    fonction_newton_Pi = (((8.d0*Pression_inc - Pression_B)/(Pression_inc + 6.d0*Pression_B)) * Pression_inc) - Pression_refl 
    RETURN
END FUNCTION fonction_newton_Pi





FUNCTION deriveFonction_newton_Pi(Pression_B , Pression_inc )
    USE Precisions
    IMPLICIT NONE
    REAL(KIND = DP) :: Pression_inc , Pression_B
    REAL(KIND = DP) :: deriveFonction_newton_Pi
    deriveFonction_newton_Pi = (((8.d0) * Pression_inc**2.d0) + 96.d0 * Pression_inc * Pression_B - 6.d0 * (Pression_B**2))&
    /((Pression_inc + 6.d0 * Pression_B)**2)
    RETURN
END FUNCTION deriveFonction_newton_Pi







FUNCTION newton_raphson_Pi(epsilon , x0_Pi , Pression_B , Pression_refl)
    USE Precisions
    IMPLICIT NONE
    REAL (KIND = DP) :: epsilon , x0_Pi , Pression_B , Pression_refl
    REAL (KIND = DP) :: newton_raphson_Pi , x1_Pi , val
    REAL (KIND = DP) , EXTERNAL :: fonction_newton_Pi , deriveFonction_newton_Pi
    Integer :: i

    val = fonction_newton_Pi(x0_Pi , Pression_B , Pression_refl)

    DO WHILE ((DABS(val) > epsilon) .OR. (i < 200))  ! on s'arrête si on a atteint la précision ou si on a fait 200 itérations (pas de convergence)
        val = fonction_newton_Pi(x0_Pi , Pression_B , Pression_refl)
        x1_Pi = x0_Pi - (fonction_newton_Pi(x0_Pi , Pression_B , Pression_refl) / deriveFonction_newton_Pi(Pression_B , x0_Pi))
        x0_Pi = x1_Pi
        i = i + 1
    END DO

    newton_raphson_Pi = x1_Pi
    RETURN
END FUNCTION newton_raphson_Pi
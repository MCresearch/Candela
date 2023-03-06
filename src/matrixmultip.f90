SUBROUTINE MULTIPAHB(M,N,K,A,LDA,B,LDB,C,LDC)
    implicit none
    INTEGER:: M,N,K,LDA,LDB,LDC
    INTEGER:: I,J,L
    COMPLEX*16:: A(LDA,*),B(LDB,*),C(LDC,*)
    COMPLEX*16:: TEMP

    DO J = 1,N
        DO I = 1,M
            TEMP = (0.0,0.0)
            DO L = 1,K
                TEMP = TEMP + DCONJG(A(L,I)) * B(L,J)
            ENDDO
            C(I,J) = TEMP
        ENDDO
    ENDDO
END SUBROUTINE
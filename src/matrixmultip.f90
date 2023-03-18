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

!down triangular matrix
SUBROUTINE DTRIMULTIPAHB(M,N,K,A,LDA,B,LDB,C,BIAS)
    implicit none
    INTEGER:: M,N,K,LDA,LDB,BIAS
    INTEGER:: I,J,L,Q
    COMPLEX*16:: A(LDA,*),B(LDB,*),C(*)
    COMPLEX*16:: TEMP

#ifdef _OPENMP
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(I,J,L,TEMP,Q)
    DO J = 1,N
        DO I = J + BIAS,M
            TEMP = (0.0,0.0)
            DO L = 1,K
                TEMP = TEMP + DCONJG(A(L,I)) * B(L,J)
            ENDDO
            Q = (J-1)*(M-BIAS-1) - (J-2)*(J-1)/2 + I-BIAS
            C(Q) = TEMP
        ENDDO
    ENDDO
    !$OMP END PARALLEL DO
#else
    Q = 1
    DO J = 1,N
        DO I = J + BIAS,M
            TEMP = (0.0,0.0)
            DO L = 1,K
                TEMP = TEMP + DCONJG(A(L,I)) * B(L,J)
            ENDDO
            C(Q) = TEMP
            Q = Q + 1;
        ENDDO
    ENDDO
#endif
END SUBROUTINE
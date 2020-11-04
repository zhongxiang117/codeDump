      BLOCK DATA BLKDATA
      COMMON /NOTHING/ TESTING
      END


      SUBROUTINE Syntax_good
      COMMON /NOTHING/ TESTING       
      END


      SUBROUTINE Culprit
      COMMON /TYPE/ Error_in_here_mistakely_identified_TYPE_as_keyword

C From here, Syntax Error: no highlights at all
      END


      SUBROUTINE Syntax_bad
      COMMON /NOTHING/ TESTING
      END


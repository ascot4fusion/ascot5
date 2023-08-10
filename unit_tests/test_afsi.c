#include "../ascot5.h"
#include <stdio.h>
#include "../boschhale.h"

int main(int argc, char** argv) {
    for(int i = 1; i <= 10000; i++) {
        printf("%le ", 1e3*i);

        for(int j = 1; j <= 4; j++) {
            printf("%le ", boschhale_sigma(j, i));
            printf("%le ", boschhale_sigmav(j, i));
        }

        printf("\n");
    }

    return 0;
}

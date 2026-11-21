/** Filename: main.cpp **/

#define STM32G431xx
#include <cstdint>
#include "stm32g4xx.h"

int main() {
    // Enable GPIOA clock
    RCC->AHB2ENR |= RCC_AHB2ENR_GPIOAEN;
    (void)RCC->AHB2ENR;

    // Configures PA12 as general purpose output
    GPIOA->MODER = (GPIOA->MODER & ~GPIO_MODER_MODE12_Msk) | (1U << GPIO_MODER_MODE12_Pos);

    // Ensures output type and low speed
    GPIOA->OTYPER  &= ~GPIO_OTYPER_OT12;
    GPIOA->OSPEEDR &= ~GPIO_OSPEEDR_OSPEED12_Msk;

    while (true) {
        // Toggle PA12 using the atomic Bit Set/Reset Register
        if (GPIOA->ODR & GPIO_ODR_OD12) {
            GPIOA->BSRR = GPIO_BSRR_BR12;
        } else {
            GPIOA->BSRR = GPIO_BSRR_BS12;
        }

        for (volatile uint32_t i = 0; i < 1'500'000; ++i) {}
    }
}

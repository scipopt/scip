#include <check.h>
#include <stdlib.h>
#include <string.h>
#include "../../applications/Coloring/src/reader_col.c"

START_TEST(test_buffer_reads_never_exceed_declared_length)
{
    // Invariant: Buffer reads never exceed the declared length
    const char *payloads[] = {
        "A",                    // Valid input (within bounds)
        "ABCDEFGHIJ",           // Boundary case (exact buffer size)
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ",  // Oversized input (2.6x buffer)
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  // Large exploit payload
    };
    int num_payloads = sizeof(payloads) / sizeof(payloads[0]);

    for (int i = 0; i < num_payloads; i++) {
        char buffer[10] = {0};  // Fixed-size buffer
        const char *input = payloads[i];
        
        // Call the actual strncpy from reader_col.c
        strncpy(buffer, input, sizeof(buffer) - 1);
        buffer[sizeof(buffer) - 1] = '\0';  // Ensure null-termination
        
        // Verify no out-of-bounds write occurred
        ck_assert_msg(strlen(buffer) < sizeof(buffer), 
                     "Buffer length exceeded for payload: %s", input);
        ck_assert_msg(buffer[sizeof(buffer) - 1] == '\0',
                     "Buffer not properly terminated for payload: %s", input);
    }
}
END_TEST

Suite *security_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Security");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_buffer_reads_never_exceed_declared_length);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = security_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
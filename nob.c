#define NOB_IMPLEMENTATION
#include "nob.h"

int main(int argc, char **argv)
{
    NOB_GO_REBUILD_URSELF(argc, argv);

    const char *program = nob_shift_args(&argc, &argv);

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "cc");
    nob_cmd_append(&cmd, "-lz", "-lpthread", "-O3");
    nob_cmd_append(&cmd, "-o", "nanodup");
    nob_cmd_append(&cmd, "nanodup.c", "thpool.c");
    if (!nob_cmd_run_sync(cmd)) return 1;

    cmd.count = 0;
    nob_cmd_append(&cmd, "cc");
    nob_cmd_append(&cmd, "-lz", "-lpthread", "-O3");
    nob_cmd_append(&cmd, "-o", "nanotrim");
    nob_cmd_append(&cmd, "nanotrim.c", "thpool.c");
    if (!nob_cmd_run_sync(cmd)) return 1;

    cmd.count = 0;
    nob_cmd_append(&cmd, "cc");
    nob_cmd_append(&cmd, "-lz", "-lpthread", "-lm", "-O3");
    nob_cmd_append(&cmd, "-o", "nanomux");
    nob_cmd_append(&cmd, "nanomux.c", "thpool.c");
    if (!nob_cmd_run_sync(cmd)) return 1;
    return 0;
}

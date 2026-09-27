// Test helper (header only; not a test of its own). Since build_v3 a malformed TESSERACT_*
// value is a hard error: the reader prints "error: NAME='value' is not valid: ..." and the
// process exits with status 2. That cannot be observed in-process, so the body runs in a
// forked child whose stderr is captured.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

// True when `body` makes the process exit with status 2 and prints the error line naming
// `flag`. Pass flag = nullptr when the body redirects stderr itself.
template <class Body>
bool exitsWithEnvError(const char* flag, Body body) {
    std::fflush(stdout);
    std::fflush(stderr);
    std::FILE* capture = std::tmpfile();
    if (!capture) return false;
    const pid_t pid = ::fork();
    if (pid < 0) { std::fclose(capture); return false; }
    if (pid == 0) {
        ::dup2(::fileno(capture), 2);
        body();
        std::fflush(stderr);
        ::_exit(0);   // the body returned: the value was accepted
    }
    int status = 0;
    ::waitpid(pid, &status, 0);
    std::string log;
    std::rewind(capture);
    char buf[4096];
    size_t n = 0;
    while ((n = std::fread(buf, 1, sizeof buf, capture)) > 0) log.append(buf, n);
    std::fclose(capture);
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 2) return false;
    return !flag || log.find(std::string("error: ") + flag + "='") != std::string::npos;
}

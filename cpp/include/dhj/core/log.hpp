// log.hpp - tee stdout and stderr into <out_dir>/run.log for the lifetime of
// the run, so the log file holds exactly what the terminal showed (the Python
// package does the same with RunLogger).
//
// Implementation: file descriptors 1 and 2 are redirected into a pipe; a
// reader thread copies everything from the pipe to the original stdout and to
// the log file. This catches every printf/fprintf in the solver without
// touching the call sites. POSIX only.
#pragma once

#include <cstdio>
#include <cstring>
#include <string>
#include <thread>

#include <fcntl.h>
#include <unistd.h>

#include "timing.hpp"

namespace dhj {

class RunLogger {
public:
    explicit RunLogger(const std::string& path, int argc, char** argv) {
        std::fflush(stdout);
        std::fflush(stderr);
        file_ = std::fopen(path.c_str(), "a");
        if (!file_) {
            std::fprintf(stderr, "warning: cannot open %s for logging; continuing without a log file\n",
                         path.c_str());
            return;
        }
        saved_out_ = ::dup(1);
        saved_err_ = ::dup(2);
        if (::pipe(fds_) != 0) {
            std::fprintf(stderr, "warning: pipe() failed; continuing without a log file\n");
            std::fclose(file_);
            file_ = nullptr;
            return;
        }
        ::dup2(fds_[1], 1);
        ::dup2(fds_[1], 2);
        ::close(fds_[1]);
        // stdout is now a pipe and would be fully buffered; keep it line-buffered
        // so progress appears promptly (and survives a crash) like a terminal.
        std::setvbuf(stdout, nullptr, _IOLBF, 0);
        active_ = true;
        reader_ = std::thread([this] { pump(); });

        std::string cmd;
        for (int i = 0; i < argc; ++i) { if (i) cmd += ' '; cmd += argv[i]; }
        char cwd[4096];
        const char* cwd_s = ::getcwd(cwd, sizeof cwd) ? cwd : "?";
        std::printf("\n%s\n", std::string(70, '#').c_str());
        std::printf("# run started %s\n# command: %s\n# cwd: %s\n", timestamp_now().c_str(), cmd.c_str(), cwd_s);
        std::printf("%s\n", std::string(70, '#').c_str());
    }

    ~RunLogger() {
        if (!active_) { if (file_) std::fclose(file_); return; }
        std::printf("# run finished %s\n", timestamp_now().c_str());
        std::fflush(stdout);
        std::fflush(stderr);
        ::dup2(saved_out_, 1);   // closes the pipe's last write ends -> reader sees EOF
        ::dup2(saved_err_, 2);
        reader_.join();
        ::close(fds_[0]);
        ::close(saved_out_);
        ::close(saved_err_);
        std::fclose(file_);
    }

    RunLogger(const RunLogger&) = delete;
    RunLogger& operator=(const RunLogger&) = delete;

private:
    void pump() {
        char buf[1 << 16];
        for (;;) {
            const ssize_t n = ::read(fds_[0], buf, sizeof buf);
            if (n <= 0) break;
            ssize_t off = 0;
            while (off < n) {
                const ssize_t w = ::write(saved_out_, buf + off, static_cast<std::size_t>(n - off));
                if (w <= 0) break;
                off += w;
            }
            std::fwrite(buf, 1, static_cast<std::size_t>(n), file_);
            std::fflush(file_);
        }
    }

    std::FILE* file_ = nullptr;
    int saved_out_ = -1, saved_err_ = -1;
    int fds_[2] = {-1, -1};
    std::thread reader_;
    bool active_ = false;
};

}  // namespace dhj

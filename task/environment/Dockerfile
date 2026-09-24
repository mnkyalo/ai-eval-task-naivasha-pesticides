# Agent-facing image.
FROM python:3.11-slim-bookworm

# --only-binary=:all: makes a missing wheel fail the build immediately with a
# clear message, instead of silently falling back to compiling scipy from
# source -- which takes far longer than any build timeout and surfaces as an
# environment-start timeout rather than as the dependency problem it is.
RUN pip install --no-cache-dir --only-binary=:all: \
    numpy==1.26.4 \
    scipy==1.11.4 \
    pandas==2.1.4

WORKDIR /app
COPY data /app/data

# Harbor starts this container and then execs into it, so it has to stay up.
# The base image's default CMD is `python3`, which reads EOF on a non-TTY
# stdin and exits at once; the container dies before Harbor can attach and
# the run fails with EnvironmentStartTimeoutError.
CMD ["sleep", "infinity"]

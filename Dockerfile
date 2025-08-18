FROM --platform=$BUILDPLATFORM rust:1-bookworm AS toolchain
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential \
      gcc-aarch64-linux-gnu g++-aarch64-linux-gnu \
      pkg-config cmake ca-certificates xxd \
 && rm -rf /var/lib/apt/lists/*

FROM toolchain AS builder
ARG TARGETARCH
COPY . .
RUN set -eux; \
    case "$TARGETARCH" in \
      amd64) RUST_TARGET=x86_64-unknown-linux-gnu; CC=gcc;  CXX=g++;  \
             export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER="$CC" ;; \
      arm64) RUST_TARGET=aarch64-unknown-linux-gnu; CC=aarch64-linux-gnu-gcc; CXX=aarch64-linux-gnu-g++; \
             export CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_LINKER="$CC" ;; \
      *) echo "Unsupported TARGETARCH: $TARGETARCH"; exit 1 ;; \
    esac; \
    rustup target add "$RUST_TARGET"; \
    export CC CXX PKG_CONFIG_ALLOW_CROSS=1; \
    export RUSTFLAGS='-C target-feature=+crt-static'; \
    export FOLDCOMP_STD_LINK=static; \
    export FOLDCOMP_LIB_DIR="$(dirname "$($CXX -print-file-name=libstdc++.a)")"; \
    cargo build --release --features foldcomp --target "$RUST_TARGET"; \
    install -D "target/$RUST_TARGET/release/folddisco" /out/folddisco

FROM debian:bookworm-slim AS runtime
COPY --from=builder /out/folddisco /usr/local/bin/folddisco
ENTRYPOINT ["/usr/local/bin/folddisco"]


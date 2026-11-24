FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y \
    build-essential gfortran \
    unzip git wget \
    bc tcl tcl-dev tcllib procps autoconf pkg-config \
    libssl-dev libcurl4-gnutls-dev libexpat1-dev \
    libxml2-dev

WORKDIR /tmp

# Install GDAL 3.12.1
RUN wget https://github.com/snale-net/pagure/archive/refs/heads/feat/updates_hdf5_netcdf.zip && \
    unzip updates_hdf5_netcdf.zip && \
    cd pagure-feat-updates_hdf5_netcdf && \
    ./pagure.sh --prefix=/build --system=ubuntu --filter=GDAL --mode=docker && \
    ./pagure.sh --prefix=/build --system=ubuntu --filter=CGAL --mode=docker

FROM fretif/python3.12-slim-spatialetl0.2.0:latest AS prebuilt
FROM python:3.12-slim-trixie AS uvbin
# Use the official Python 3.12 slim image for a lightweight base.
FROM python:3.12-slim AS python-builder

# Prevent Python from writing .pyc files to disk and enable unbuffered logging.
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install build dependencies (e.g., gcc) needed to compile Python packages.
# Note: We remove package lists afterwards to keep the image small.
RUN apt-get update && apt-get install --no-install-recommends -y \
    gcc g++ pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Copy the compiled assets that we need to run from the builder image
# System
COPY --from=prebuilt /lib/x86_64-linux-gnu/libnghttp2.so.14 /usr/local/lib/libnghttp2.so.14
COPY --from=prebuilt /lib/x86_64-linux-gnu/libnghttp3.so.9 /usr/local/lib/libnghttp3.so.9
COPY --from=prebuilt /lib/x86_64-linux-gnu/libngtcp2_crypto_gnutls.so.8 /usr/local/lib/libngtcp2_crypto_gnutls.so.8
COPY --from=prebuilt /lib/x86_64-linux-gnu/libngtcp2.so.16 /usr/local/lib/libngtcp2.so.16
COPY --from=prebuilt /lib/x86_64-linux-gnu/libidn2.so.0 /usr/local/lib/libidn2.so.0
COPY --from=prebuilt /lib/x86_64-linux-gnu/libgnutls.so.30 /usr/local/lib/libgnutls.so.30
COPY --from=prebuilt /lib/x86_64-linux-gnu/libgssapi_krb5.so.2 /usr/local/lib/libgssapi_krb5.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/libunistring.so.5 /usr/local/lib/libunistring.so.5
COPY --from=prebuilt /lib/x86_64-linux-gnu/librtmp.so.1 /usr/local/lib/librtmp.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libssh2.so.1 /usr/local/lib/libssh2.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libpsl.so.5 /usr/local/lib/libpsl.so.5
COPY --from=prebuilt /lib/x86_64-linux-gnu/libldap.so.2 /usr/local/lib/libldap.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/liblber.so.2 /usr/local/lib/liblber.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/libbrotlidec.so.1 /usr/local/lib/libbrotlidec.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libexpat.so.1 /usr/local/lib/libexpat.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libgfortran.so.5 /usr/local/lib/libgfortran.so.5
COPY --from=prebuilt /lib/x86_64-linux-gnu/libcurl-gnutls.so.4 /usr/local/lib/libcurl-gnutls.so.4
COPY --from=prebuilt /lib/x86_64-linux-gnu/libsasl2.so.2 /usr/local/lib/libsasl2.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/libbrotlicommon.so.1 /usr/local/lib/libbrotlicommon.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libp11-kit.so.0 /usr/local/lib/libp11-kit.so.0
COPY --from=prebuilt /lib/x86_64-linux-gnu/libtasn1.so.6 /usr/local/lib/libtasn1.so.6
COPY --from=prebuilt /lib/x86_64-linux-gnu/libkrb5.so.3 /usr/local/lib/libkrb5.so.3
COPY --from=prebuilt /lib/x86_64-linux-gnu/libk5crypto.so.3 /usr/local/lib/libk5crypto.so.3
COPY --from=prebuilt /lib/x86_64-linux-gnu/libcom_err.so.2 /usr/local/lib/libcom_err.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/libkrb5support.so.0 /usr/local/lib/libkrb5support.so.0
COPY --from=prebuilt /lib/x86_64-linux-gnu/libkeyutils.so.1 /usr/local/lib/libkeyutils.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libxml2.so.2 /usr/local/lib/libxml2.so.2
#COPY --from=prebuilt /lib/x86_64-linux-gnu/libwebp.so.7 /usr/local/lib/libwebp.so.7
#COPY --from=prebuilt /lib/x86_64-linux-gnu/libjpeg.so.62 /usr/local/lib/libjpeg.so.62
#COPY --from=prebuilt /lib/x86_64-linux-gnu/libdeflate.so.0 /usr/local/lib/libdeflate.so.0
#COPY --from=prebuilt /lib/x86_64-linux-gnu/libLerc.so.4 /usr/local/lib/libLerc.so.4
#COPY --from=prebuilt /lib/x86_64-linux-gnu/libjbig.so.0 /usr/local/lib/libjbig.so.0

# Lib
COPY --from=prebuilt /build/zlib/gcc142/1.2.11/lib /usr/local/lib
COPY --from=prebuilt /build/lapack-blas/gcc142/3.9.1/lib /usr/local/lib
COPY --from=prebuilt /build/hdf5/gcc142/1.14.6/lib /usr/local/lib
COPY --from=prebuilt /build/netcdf/hdf5.146/gcc142/c/4.9.3/lib /usr/local/lib
COPY --from=prebuilt /build/tiff/gcc142/4.4.0/lib /usr/local/lib
COPY --from=prebuilt /build/geos/gcc142/3.14.1/lib /usr/local/lib
COPY --from=prebuilt /build/sqlite/gcc142/3.36.0/lib /usr/local/lib
COPY --from=prebuilt /build/proj/gcc142/9.7.1/lib /usr/local/lib
COPY --from=prebuilt /build/gdal/gcc142/3.12.1/lib /usr/local/lib
COPY --from=prebuilt /build/boost/py37/gcc142/1.90.0/lib /usr/local/lib
COPY --from=prebuilt /build/cgal/gcc142/6.1/lib /usr/local/lib
RUN ldconfig

# Include
COPY --from=prebuilt /build/gdal/gcc142/3.12.1/include /usr/local/include
COPY --from=prebuilt /build/boost/py37/gcc142/1.90.0/include /usr/local/include
COPY --from=prebuilt /build/cgal/gcc142/6.1/include /usr/local/include

# Bin
COPY --from=prebuilt /build/gdal/gcc142/3.12.1/bin/ /usr/local/bin/

ENV CGAL_DIR=/usr/local/include
ENV UV_NO_DEV=1
ENV UV_SYSTEM_PYTHON=1
ENV UV_COMPILE_BYTECODE=1

# Set the working directory to /app
WORKDIR /app

# Copy the entire project into the container.
COPY . .

RUN --mount=from=ghcr.io/astral-sh/uv,source=/uv,target=/bin/uv \
 uv sync --extra interpcgal --extra gdal --extra symphonie --locked --no-editable

FROM python:3.12-slim AS prod

# Copy the compiled assets that we need to run from the builder image
# System
COPY --from=prebuilt /lib/x86_64-linux-gnu/libnghttp2.so.14 /usr/local/lib/libnghttp2.so.14
COPY --from=prebuilt /lib/x86_64-linux-gnu/libnghttp3.so.9 /usr/local/lib/libnghttp3.so.9
COPY --from=prebuilt /lib/x86_64-linux-gnu/libngtcp2_crypto_gnutls.so.8 /usr/local/lib/libngtcp2_crypto_gnutls.so.8
COPY --from=prebuilt /lib/x86_64-linux-gnu/libngtcp2.so.16 /usr/local/lib/libngtcp2.so.16
COPY --from=prebuilt /lib/x86_64-linux-gnu/libidn2.so.0 /usr/local/lib/libidn2.so.0
COPY --from=prebuilt /lib/x86_64-linux-gnu/libgnutls.so.30 /usr/local/lib/libgnutls.so.30
COPY --from=prebuilt /lib/x86_64-linux-gnu/libgssapi_krb5.so.2 /usr/local/lib/libgssapi_krb5.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/libunistring.so.5 /usr/local/lib/libunistring.so.5
COPY --from=prebuilt /lib/x86_64-linux-gnu/librtmp.so.1 /usr/local/lib/librtmp.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libssh2.so.1 /usr/local/lib/libssh2.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libpsl.so.5 /usr/local/lib/libpsl.so.5
COPY --from=prebuilt /lib/x86_64-linux-gnu/libldap.so.2 /usr/local/lib/libldap.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/liblber.so.2 /usr/local/lib/liblber.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/libbrotlidec.so.1 /usr/local/lib/libbrotlidec.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libexpat.so.1 /usr/local/lib/libexpat.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libgfortran.so.5 /usr/local/lib/libgfortran.so.5
COPY --from=prebuilt /lib/x86_64-linux-gnu/libcurl-gnutls.so.4 /usr/local/lib/libcurl-gnutls.so.4
COPY --from=prebuilt /lib/x86_64-linux-gnu/libsasl2.so.2 /usr/local/lib/libsasl2.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/libbrotlicommon.so.1 /usr/local/lib/libbrotlicommon.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libp11-kit.so.0 /usr/local/lib/libp11-kit.so.0
COPY --from=prebuilt /lib/x86_64-linux-gnu/libtasn1.so.6 /usr/local/lib/libtasn1.so.6
COPY --from=prebuilt /lib/x86_64-linux-gnu/libkrb5.so.3 /usr/local/lib/libkrb5.so.3
COPY --from=prebuilt /lib/x86_64-linux-gnu/libk5crypto.so.3 /usr/local/lib/libk5crypto.so.3
COPY --from=prebuilt /lib/x86_64-linux-gnu/libcom_err.so.2 /usr/local/lib/libcom_err.so.2
COPY --from=prebuilt /lib/x86_64-linux-gnu/libkrb5support.so.0 /usr/local/lib/libkrb5support.so.0
COPY --from=prebuilt /lib/x86_64-linux-gnu/libkeyutils.so.1 /usr/local/lib/libkeyutils.so.1
COPY --from=prebuilt /lib/x86_64-linux-gnu/libxml2.so.2 /usr/local/lib/libxml2.so.2

# Lib
COPY --from=prebuilt /build/zlib/gcc142/1.2.11/lib /usr/local/lib
COPY --from=prebuilt /build/lapack-blas/gcc142/3.9.1/lib /usr/local/lib
COPY --from=prebuilt /build/hdf5/gcc142/1.14.6/lib /usr/local/lib
COPY --from=prebuilt /build/netcdf/hdf5.146/gcc142/c/4.9.3/lib /usr/local/lib
COPY --from=prebuilt /build/tiff/gcc142/4.4.0/lib /usr/local/lib
COPY --from=prebuilt /build/geos/gcc142/3.14.1/lib /usr/local/lib
COPY --from=prebuilt /build/sqlite/gcc142/3.36.0/lib /usr/local/lib
COPY --from=prebuilt /build/proj/gcc142/9.7.1/lib /usr/local/lib
COPY --from=prebuilt /build/gdal/gcc142/3.12.1/lib /usr/local/lib
COPY --from=prebuilt /build/boost/py37/gcc142/1.90.0/lib /usr/local/lib
COPY --from=prebuilt /build/cgal/gcc142/6.1/lib /usr/local/lib
RUN ldconfig

# Share
COPY --from=prebuilt /build/proj/gcc142/9.7.1/share/ /usr/local/share/
COPY --from=prebuilt /build/gdal/gcc142/3.12.1/share/ /usr/local/share/

# Bin
COPY --from=prebuilt /build/gdal/gcc142/3.12.1/bin/ /usr/local/bin/

# Export
ENV GDAL_DATA="/usr/local/share/gdal"

# Copy the environment, but not the source code
COPY --from=python-builder --chown=app:app /app/.venv /app/.venv

ENV PATH="/app/.venv/bin:$PATH"


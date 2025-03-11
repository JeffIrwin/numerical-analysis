
FROM rockylinux:9

WORKDIR /workdir

RUN dnf update -y
RUN dnf install -y git
#RUN dnf install -y gfortran  # gfortran 11

RUN dnf install -y gcc-toolset-13-gcc-gfortran
ENV PATH="$PATH:/opt/rh/gcc-toolset-13/root/usr/bin/"

RUN gfortran --version

ADD https://github.com/fortran-lang/fpm/releases/download/current/fpm-linux-x86_64 ./fpm
RUN chmod +x fpm
RUN mv fpm /usr/local/bin
RUN fpm --version

ARG BRANCH="main"
RUN echo "BRANCH = $BRANCH"

RUN git clone https://github.com/jeffirwin/numerical-analysis --branch "$BRANCH"
WORKDIR /workdir/numerical-analysis

RUN fpm test --profile debug
RUN fpm test --profile release


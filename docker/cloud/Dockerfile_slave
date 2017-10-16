from ubuntu:16.04

ARG DOCKER_ENGINE_VERSION
ARG DOCKER_COMPOSE_VERSION
ARG JENKINS_SSH_PUBLIC_KEY

ENV DOCKER_ENGINE_VERSION=${DOCKER_ENGINE_VERSION} \
    DOCKER_COMPOSE_VERSION=${DOCKER_COMPOSE_VERSION}

# Install SSH server, Git, and cURL
RUN apt-get update && apt-get install -y \
      openssh-server \
      git \
      curl \
      default-jdk \
      && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir -p /var/run/sshd

# Get Docker Engine binaries
RUN curl -fsSLO https://get.docker.com/builds/Linux/x86_64/docker-${DOCKER_ENGINE_VERSION}.tgz && \
    tar --strip-components=1 -xvzf docker-${DOCKER_ENGINE_VERSION}.tgz -C /usr/local/bin && \
    rm docker-${DOCKER_ENGINE_VERSION}.tgz && \
    docker --version

# Install Docker Compose
RUN curl -fsSL https://github.com/docker/compose/releases/download/${DOCKER_COMPOSE_VERSION}/docker-compose-`uname -s`-`uname -m` > /usr/local/bin/docker-compose && \
    chmod +x /usr/local/bin/docker-compose && \
    docker-compose --version

# Copy SSH Jenkins credentials
COPY ${JENKINS_SSH_PUBLIC_KEY} /root/.ssh/authorized_keys

EXPOSE 22
VOLUME /var/run/docker.sock

CMD ["/usr/sbin/sshd", "-D"]

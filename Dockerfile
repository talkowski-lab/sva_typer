FROM rust:1.86

# WORKDIR /usr/src
#
# RUN git clone https://github.com/talkowski-lab/sva_typer.git 
#
# RUN cd sva_typer && cargo install --path .

WORKDIR /usr/src/sva_typer

COPY . .

RUN cargo install --path .

CMD ["sva_typer"]

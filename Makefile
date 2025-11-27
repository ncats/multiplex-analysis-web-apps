# Usage:
# make insert_image_metadata ENV_NAME=leandro-robert ENV_PY_VER=3.12 DATE=2025-11-26 BUILD_VER=01

insert_metadata:
	@IMAGE_NAME=mawa-frontend; \
	TAG="$(DATE)-v$(BUILD_VER)-$(ENV_NAME)-env"; \
	ENV_FILE="environment-$(ENV_NAME)-compatible.yml"; \
	WHO_ADDED="andrewweisman"; \
	DIGEST=$$(docker inspect --format='{{.Id}}' frontend:$$TAG); \
	COMMIT=$$(git rev-parse HEAD); \
	SQL="insert into general_schema.image_metadata_table (image_id, name, tag, git_commit, environment_yaml_file, archive_compatibility_id, who_added) values ('$$DIGEST', '$$IMAGE_NAME', '$$TAG', '$$COMMIT', '$$ENV_FILE', '$(ENV_NAME)', '$$WHO_ADDED');"; \
	echo "$$SQL"; \
    docker exec -i postgres psql "postgresql://postgres:password@database:5432/app_a_app_db" -c "$$SQL"

build:
	@ENV_NAME=$(ENV_NAME) ENV_PY_VER=$(ENV_PY_VER) DATE=$(DATE) BUILD_VER=$(BUILD_VER) docker compose build

up:
	@ENV_NAME=$(ENV_NAME) ENV_PY_VER=$(ENV_PY_VER) DATE=$(DATE) BUILD_VER=$(BUILD_VER) docker compose up

build_up:
	@ENV_NAME=$(ENV_NAME) ENV_PY_VER=$(ENV_PY_VER) DATE=$(DATE) BUILD_VER=$(BUILD_VER) docker compose up --build

down:
	@ENV_NAME=$(ENV_NAME) ENV_PY_VER=$(ENV_PY_VER) DATE=$(DATE) BUILD_VER=$(BUILD_VER) docker compose down

build_insert_up:
	@echo "Building images..."
	@ENV_NAME=$(ENV_NAME) ENV_PY_VER=$(ENV_PY_VER) DATE=$(DATE) BUILD_VER=$(BUILD_VER) docker compose build
	@echo "Starting postgres..."
	@ENV_NAME=$(ENV_NAME) ENV_PY_VER=$(ENV_PY_VER) DATE=$(DATE) BUILD_VER=$(BUILD_VER) docker compose up -d postgres
	@sleep 3
	@echo "Inserting metadata..."
	@$(MAKE) insert_metadata ENV_NAME=$(ENV_NAME) DATE=$(DATE) BUILD_VER=$(BUILD_VER)
	@echo "Stopping postgres..."
	@ENV_NAME=$(ENV_NAME) ENV_PY_VER=$(ENV_PY_VER) DATE=$(DATE) BUILD_VER=$(BUILD_VER) docker compose down
	@echo "Starting all services..."
	@ENV_NAME=$(ENV_NAME) ENV_PY_VER=$(ENV_PY_VER) DATE=$(DATE) BUILD_VER=$(BUILD_VER) docker compose up

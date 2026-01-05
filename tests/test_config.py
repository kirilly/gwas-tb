"""Tests for configuration module."""

import pytest

from src.config import Config


def test_config_from_yaml(temp_config_file):
    """Test loading config from YAML."""
    config = Config.from_yaml(temp_config_file)

    assert config.project.name == "test-gwas"
    assert config.project.random_seed == 42
    assert config.qc.min_maf == 0.01


def test_config_defaults():
    """Test default config values."""
    config = Config()

    assert config.project.name == "tb-gwas"
    assert config.project.random_seed == 42
    assert config.qc.min_maf == 0.01
    assert config.phylogeny.method == "fasttree"


def test_config_validation():
    """Test config validation."""
    config = Config()
    errors = config.validate()
    assert len(errors) == 0

    # Invalid MAF
    config.qc.min_maf = -0.1
    errors = config.validate()
    assert len(errors) > 0
    assert "MAF" in errors[0]


def test_config_to_dict():
    """Test config serialization."""
    config = Config()
    d = config.to_dict()

    assert "project" in d
    assert d["project"]["name"] == "tb-gwas"


def test_config_file_not_found():
    """Test error on missing config file."""
    with pytest.raises(FileNotFoundError):
        Config.from_yaml("/nonexistent/config.yaml")

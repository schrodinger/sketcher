## Schrödinger 2D Sketcher

![build-status](https://github.com/schrodinger/sketcher/actions/workflows/nightly-build.yml/badge.svg?branch=main)

[![sketcher](https://github.com/schrodinger/sketcher/blob/main/.github/schrodinger-sketcher-screenshot.png)](https://www.schrodinger.com/2dsketcher)

The Schrödinger 2D Sketcher is an open-source application for drawing and editing chemical structures. The entire sketcher is built on top of the RDKit open source toolkit, which serves as its underlying chemical model. It provides a user-friendly interface for creating molecules, reactions, and other chemical diagrams, which can be used independently or integrated into other cheminformatics workflows.

This project is released by Schrödinger, Inc. and is available under an open-source license to foster collaboration and development within the scientific community.

## Features

* Intuitive drawing tools for atoms, bonds, rings, and functional groups.
* Support for many chemical file formats.
* Cleanup and layout algorithms for generating clear 2D representations.
* Support for stereochemistry representation.
* Integration capabilities for use in other applications, including as a web component.

## Access and Training

**[Open Access Online Version](https://www.schrodinger.com/2dsketcher)** -- Explore the sketcher directly in your web browser.

**[Training and Video Walkthrough](https://www.schrodinger.com/sites/default/files/s3/public/2D-Sketcher/2023-2/Content/Resources/Videos/2D_Sketcher.mp4)** -- Learn how to use the sketcher with this guided video.

## Build Prerequisites

* A C++ compiler supporting C++20 or later
* CMake (version 3.24 or later)
* All [required dependencies](https://github.com/schrodinger/sketcher/blob/main/external/versions.json) installed and accessible to CMake

## Support and Community

For questions, support, or discussions related to the Schrödinger 2D Sketcher, please use the [GitHub issue tracker](https://github.com/schrodinger/sketcher/issues).

To contribute code, please follow these steps:

* **Fork the repository:** Fork the `schrodinger/sketcher` repository to your own GitHub account.
* **Create a branch:** Create a new branch in your forked repository for your changes.
* **Implement your changes:** Make your code changes or additions.
* **Test your changes:** Ensure your changes build correctly and pass any existing tests. Add new tests for new features if applicable.
* **Submit a Pull Request:** Open a pull request from your branch to the `main` branch of the `schrodinger/sketcher` repository. Please provide a clear description of your changes.

All contributions are subject to the terms of the BSD license.

## License

The code is released under the [BSD 3-Clause License](https://github.com/schrodinger/sketcher/blob/master/LICENSE).

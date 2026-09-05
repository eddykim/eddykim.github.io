---
title: Starting This Blog
lang: en
lang-exclusive: ["en"]
permalink: /posts/start-blog/
page_id: start-blog
date: 2026-09-05 13:00:00 +0900
categories: [Notes]
tags: [intro, blog, jekyll, github-pages]
description: "Why I am starting a research notebook on optical metrology, and what broke while setting it up on GitHub Pages."
---

## Why start a blog

After finishing a long PhD, I looked back and found that everything I had studied, researched, and coded was scattered around in no particular order. Going through it again, my honest reaction was that I had studied less than I thought, and what I did study I never really finished properly. For eight years I built optical metrology instruments — spectroscopic ellipsometry, Mueller matrix polarimetry, AI-assisted metrology, spectroscopic reflectometry — from the hardware through alignment and calibration to the analysis software. Almost all of it ended up locked inside papers, patents, and lab code, or scattered across slide decks and notes.

I have always admired people who keep good technical blogs. Writing well is hard, and explaining something clearly is harder. But the material is there, and with an LLM to help I thought I might be able to turn it into something worth reading. Publishing the process itself, and writing down what I learn along the way — if I keep filling in those gaps, that alone becomes a kind of credibility in this field.

## Why GitHub Pages

I looked at Velog first. It has a low barrier to entry and a healthy developer community. In the end I went with GitHub Pages and Jekyll.

- My code and profile already live at `github.com/eddykim`, so keeping the blog on the same account and the same workflow (git commit, PR) felt natural.
- Posts are just markdown files, which means a local editor and git version control. That is the familiar path for an engineer.
- I write a lot about metrology, which means equations, code, and figures. I needed room to customize.
- This one mattered more than I expected: I use Obsidian, so drafting and editing markdown before publishing is easy. A well-written Obsidian note can become a post with only minor edits.

For the theme I picked [Chirpy](https://github.com/cotes2020/jekyll-theme-chirpy). It is actively maintained and reasonably well documented, and the `chirpy-starter` template saved me from unzipping a theme and shuffling files around by hand.

## Things that broke while setting it up

Setting up the blog turned out to be the first assignment. Here is what cost me a few days.

The one that took longest was posts not showing up. I pushed, and the site stayed exactly the same. My first thought was that I had messed up a branch, so I went digging through the git log. Nothing wrong there. The Actions log said the build succeeded, which made it more confusing. It turned out that one reference link in a post started with `http://`, and the link checker flagged it. The workflow runs build, then check, then deploy — and since it stopped at the check, deploy never even started. One link kept six posts offline for three days.

Scheduled publishing was another one. Jekyll will not publish a post dated in the future. But GitHub Pages only builds when you push. Put those two together and a scheduled post never appears, because nobody triggers a build once the time passes. It only started working after I added a job that builds on a timer.

I also found out late that a file I never meant to publish was public. A working-notes file at the top of the repository was not in the exclude list, so it was served as part of the site. Opening the URL directly returned it. Jekyll copies everything you do not explicitly exclude, and I had forgotten that.

Then there was the time I fixed something and the page refused to change. Chirpy caches pages in the browser through a service worker and shows the cached copy first, so any page you have already visited keeps serving the old version. As it happens, the theme already ships an "update available" notice — though at the time deployment was stuck anyway, so there was nothing for it to announce.

The most useful thing to come out of all this was not any single fix. It was setting things up so I can run exactly the same checks the server runs, on my own machine. Two commands before pushing, and none of the problems above make it out of my laptop.

## What I plan to write about

- **Optics** — from electromagnetic fundamentals to the Fresnel equations, then polarization and Mueller matrices, then multilayer thin-film modeling. Ellipsometry and spectroscopic reflectometry are the center of gravity.
- **Computation** — how you get a physical quantity out of measured data. Starting from least squares and gradient descent, through Gauss-Newton, Levenberg-Marquardt, and global optimization. Image processing and machine learning come next. What interests me most is how physics-based forward models and data-driven methods fit together.
- **Surface analysis** — a new thread. XPS, starting from the photoelectric effect and binding energy, then why it is surface sensitive, how to read the structure of a spectrum, and how far you can trust a quantification.
- **Research and tooling notes** — writing up ongoing work (a line-scan spectroscopic reflectometer, among others) in a form you can actually run, plus records like the one above of things that broke.

Code that appears in a post is committed in runnable form alongside it. Posts carry the essential part inline and link to the rest.

Posts are written in Korean first. English versions like this one are added for the topics likely to be useful outside Korea.

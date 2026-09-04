# frozen_string_literal: true

source "https://rubygems.org"

gem "jekyll-theme-chirpy", "~> 7.6"

# 다국어(ko/en) 빌드. 워크플로가 `bundle exec jekyll b` 로 빌드하므로 사용 가능하다.
gem "jekyll-polyglot", "~> 1.14"

gem "html-proofer", "~> 5.0", group: :test

platforms :windows, :jruby do
  gem "tzinfo", ">= 1", "< 3"
  gem "tzinfo-data"
end

gem "wdm", "~> 0.2.0", :platforms => [:windows]

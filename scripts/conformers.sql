-- phpMyAdmin SQL Dump
-- version 4.0.10deb1
-- http://www.phpmyadmin.net
--
-- Host: porky
-- Generation Time: Feb 22, 2015 at 09:19 AM
-- Server version: 5.5.41-0ubuntu0.14.04.1
-- PHP Version: 5.5.9-1ubuntu4.6

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `conformers`
--
CREATE DATABASE IF NOT EXISTS `conformers` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci;
USE `conformers`;

-- --------------------------------------------------------

--
-- Table structure for table `names`
--

CREATE TABLE IF NOT EXISTS `names` (
  `smile` varchar(256) NOT NULL,
  `name` varchar(64) NOT NULL,
  PRIMARY KEY (`smile`,`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `structures`
--

CREATE TABLE IF NOT EXISTS `structures` (
  `smile` varchar(256) NOT NULL,
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT,
  `weight` float NOT NULL,
  `nconfs` int(11) NOT NULL DEFAULT '0',
  `sdfloc` varchar(256) DEFAULT NULL COMMENT 'file path to sdf.gz file',
  PRIMARY KEY (`smile`),
  UNIQUE KEY `id` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;

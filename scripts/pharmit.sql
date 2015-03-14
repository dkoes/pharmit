-- phpMyAdmin SQL Dump
-- version 4.0.10deb1
-- http://www.phpmyadmin.net
--
-- Host: localhost
-- Generation Time: Mar 14, 2015 at 04:53 PM
-- Server version: 5.5.41-0ubuntu0.14.04.1
-- PHP Version: 5.5.9-1ubuntu4.6

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `pharmit`
--
CREATE DATABASE IF NOT EXISTS `pharmit` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci;
USE `pharmit`;

-- --------------------------------------------------------

--
-- Table structure for table `databases`
--

CREATE TABLE IF NOT EXISTS `databases` (
  `email` varchar(256) NOT NULL,
  `name` text NOT NULL,
  `description` text NOT NULL,
  `id` varchar(64) NOT NULL,
  `isprivate` tinyint(1) NOT NULL DEFAULT '0',
  `status` varchar(32) NOT NULL,
  `message` text NOT NULL,
  `directory` text NOT NULL,
  `submitted` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `completed` timestamp NOT NULL DEFAULT '0000-00-00 00:00:00',
  `nummols` int(11) NOT NULL DEFAULT '0',
  `numconfs` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  UNIQUE KEY `id` (`id`),
  KEY `email` (`email`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `users`
--

CREATE TABLE IF NOT EXISTS `users` (
  `email` varchar(256) NOT NULL,
  `password` text NOT NULL,
  `name` text NOT NULL,
  `institution` text NOT NULL,
  `maxprivatedbs` int(11) NOT NULL DEFAULT '1',
  `maxprivateconfs` int(11) NOT NULL DEFAULT '1000000',
  `maxdbs` int(11) NOT NULL DEFAULT '0',
  `maxconfs` int(11) NOT NULL DEFAULT '10000000',
  PRIMARY KEY (`email`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
